#pragma ident "$Id: ComputeStaClock.cpp $"

/**
 * @file ComputeStaClock.cpp
 * This is a class to compute the station clock bias.
 */

//============================================================================
//
//  This file is part of GPSTk, the GPS Toolkit.
//
//  The GPSTk is free software; you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published
//  by the Free Software Foundation; either version 2.1 of the License, or
//  any later version.
//
//  The GPSTk is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with GPSTk; if not, write to the Free Software Foundation,
//  Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//  Kaifa Kuang - Wuhan University. 2017
//
//============================================================================


#include "ComputeStaClock.hpp"

using namespace std;

namespace gpstk
{

    // Return a string identifying this object.
    std::string ComputeStaClock::getClassName() const
    { return "ComputeStaClock"; }


    /** Return a reference to a satTypeValueMap object after differencing
     *  data type values given in 'diffTypes' field with the previous
     *  corresponding type.
     *
     * @param gData      Data object holding the data.
     */
    gnssSatTypeValue& ComputeStaClock::Process(gnssSatTypeValue& gData)
        throw(ProcessingException)
    {
        try
        {
            // Build a gnssRinex object and fill it with data
            gnssRinex g1;
            g1.header = gData.header;
            g1.body = gData.body;

            // Call the Process() method with the appropriate input object
            Process(g1);

            // Update the original gnssSatTypeValue object with the results
            gData.body = g1.body;

            return gData;
        }
        catch(Exception& u)
        {
            // Throw an exception if something unexpected happens
            ProcessingException e( getClassName() + u.what() );
            GPSTK_THROW(e);
        }
    }

    /** Return a reference to a satTypeValueMap object after differencing
     *  data type values given in 'diffTypes' field with the previous
     *  corresponding type.
     *
     * @param gData      Data object holding the data.
     */
    gnssRinex& ComputeStaClock::Process(gnssRinex& gData)
        throw(ProcessingException)
    {
        try
        {
            Process(gData.body);
            return gData;
        }
        catch(Exception& u)
        {
            // Throw an exception if something unexpected happens
            ProcessingException e( getClassName() + u.what() );
            GPSTK_THROW(e);
        }
    }


     /** Return a gnssDataMap object, adding the new data generated
      *  when calling this object.
      *
      * @param gData    Data object holding the data.
      */
    gnssDataMap& ComputeStaClock::Process(gnssDataMap& gData)
        throw(ProcessingException)
    {
        SourceIDSet sourceRejectedSet;

        sourceClockMap.clear();

        for(gnssDataMap::iterator gdmIt = gData.begin();
            gdmIt != gData.end();
            ++gdmIt)
        {
            for(sourceDataMap::iterator sdmIt = gdmIt->second.begin();
                sdmIt != gdmIt->second.end();
                ++sdmIt)
            {
                SourceID source( sdmIt->first );

//                cout << setw(5) << source.sourceName.substr(0,4) << endl;

                sourceRejected = false;

                Process( sdmIt->second );

                if(sourceRejected)
                {
                    sourceRejectedSet.insert(source);
                    continue;
                }
                else
                {
                    sourceClockMap[source] = clock;
                }
            }
        }

        gData.removeSourceID(sourceRejectedSet);

        return gData;

    }  // End of method 'ComputeStaClock::Process()'


    /** Return a reference to a satTypeValueMap object after differencing
     *  data type values given in 'diffTypes' field with the previous
     *  corresponding type.
     *
     * @param gData      Data object holding the data.
     */
    satTypeValueMap& ComputeStaClock::Process(satTypeValueMap& gData)
        throw(ProcessingException)
    {
        try
        {
            SatIDSet gpsRejectedSet;
            SatIDSet galRejectedSet;

            Vector<SatID> allSatVec( gData.getVectorOfSatID() );

            vector<SatID> gpsSatVec;
            vector<SatID> galSatVec;

            for(int i=0; i<allSatVec.size(); ++i)
            {
                SatID sat( allSatVec[i] );

                if(sat.system == SatID::systemGPS)
                {
                    gpsSatVec.push_back( sat );
                }
                else if(sat.system == SatID::systemGalileo)
                {
                    galSatVec.push_back( sat );
                }
            }

            int total(0);

            // gps
            total = gpsSatVec.size();

            if(total != 0)
            {
                int rejected( 0 );
                int reminder( 0 );

                bool rejectAll( false );

//                cout << setw(5) << total << endl;

                while( true )
                {
                    rejected = gpsRejectedSet.size();
                    reminder = total - rejected;

                    if(reminder < minSatGPS)
                    {
                        rejectAll = true;
                        break;
                    }

                    Vector<double> prefit(total,0.0);
                    Vector<double> weight(total,0.0);

                    for(int i=0; i<total; ++i)
                    {
                        SatID sat( gpsSatVec[i] );

                        if(gpsRejectedSet.find(sat) != gpsRejectedSet.end())
                            continue;

                        prefit(i) = gData.getValue(sat,TypeID::prefitC12);
                        weight(i) = gData.getValue(sat,TypeID::weight);
                    }

                    clock = dot(weight,prefit)/sum(weight);

                    Vector<double> postfit(total,0.0);

                    for(int i=0; i<total; ++i)
                    {
                        SatID sat( gpsSatVec[i] );

                        if(gpsRejectedSet.find(sat) != gpsRejectedSet.end())
                            continue;

                        postfit(i) = prefit(i) - clock;
                    }

                    double rms( norm(postfit)/std::sqrt(reminder) );

                    if(rms <= maxRMSGPS)
                    {
                        rejectAll = false; break;
                    }
                    else
                    {
                        int pos( 0 );
                        double value( 0.0 );

                        for(int i=0; i<total; ++i)
                        {
                            if(std::abs(postfit(i)) > value)
                            {
                                pos = i;
                                value = std::abs( postfit(i) );
                            }
                        }

                        rejectAll = false;
                        gpsRejectedSet.insert( gpsSatVec[pos] );
                    }
                }

                if( rejectAll )
                {
                    sourceRejected = true;
                    return gData;
                }

                gData.removeSatID( gpsRejectedSet );
            }


            // galileo
            total = galSatVec.size();

            if(total != 0)
            {
                int rejected( 0 );
                int reminder( 0 );

                bool rejectAll( false );

//                cout << setw(5) << total << endl;

                while( true )
                {
                    rejected = galRejectedSet.size();
                    reminder = total - rejected;

                    if(reminder < minSatGal)
                    {
                        rejectAll = true;
                        break;
                    }

                    Vector<double> prefit(total,0.0);
                    Vector<double> weight(total,0.0);

                    for(int i=0; i<total; ++i)
                    {
                        SatID sat( galSatVec[i] );

                        if(galRejectedSet.find(sat) != galRejectedSet.end())
                            continue;

                        prefit(i) = gData.getValue(sat,TypeID::prefitC12);
                        weight(i) = gData.getValue(sat,TypeID::weight);
                    }

                    clock = dot(weight,prefit)/sum(weight);

                    Vector<double> postfit(total,0.0);

                    for(int i=0; i<total; ++i)
                    {
                        SatID sat( galSatVec[i] );

                        if(galRejectedSet.find(sat) != galRejectedSet.end())
                            continue;

                        postfit(i) = prefit(i) - clock;
                    }

                    double rms( norm(postfit)/std::sqrt(reminder) );

                    if(rms <= maxRMSGal)
                    {
                        rejectAll = false; break;
                    }
                    else
                    {
                        int pos( 0 );
                        double value( 0.0 );

                        for(int i=0; i<total; ++i)
                        {
                            if(std::abs(postfit(i)) > value)
                            {
                                pos = i;
                                value = std::abs( postfit(i) );
                            }
                        }

                        rejectAll = false;
                        galRejectedSet.insert( galSatVec[pos] );
                    }
                }

//                cout << "rejectAll: " << setw(5) << rejectAll << endl;
                if( rejectAll )
                {
                    sourceRejected = true;
                    return gData;
                }

                gData.removeSatID( galRejectedSet );
            }

//            allSatVec = gData.getVectorOfSatID();
//
//            gpsSatVec.clear();
//            galSatVec.clear();
//            for(int i=0; i<allSatVec.size(); ++i)
//            {
//                SatID sat( allSatVec[i] );
//
//                if(sat.system == SatID::systemGPS)
//                {
//                    gpsSatVec.push_back( sat );
//                }
//                else if(sat.system == SatID::systemGalileo)
//                {
//                    galSatVec.push_back( sat );
//                }
//            }
//
//            Vector<double> prefit;
//            Vector<double> weight;
//
//            // gps clock
//            prefit.resize(gpsSatVec.size(), 0.0);
//            weight.resize(gpsSatVec.size(), 0.0);
//
//            for(int i=0; i<gpsSatVec.size(); ++i)
//            {
//                SatID sat( gpsSatVec[i] );
//
//                prefit(i) = gData.getValue(sat,TypeID::prefitC12);
//                weight(i) = gData.getValue(sat,TypeID::weight);
//            }
//
//            clock = dot(weight,prefit)/sum(weight);
//
//            for(int i=0; i<gpsSatVec.size(); ++i)
//            {
//                SatID sat( gpsSatVec[i] );
//
//                gData[sat][TypeID::cdtSta] = clock;
//            }

            // galileo clock
//            prefit.resize(galSatVec.size(), 0.0);
//            weight.resize(galSatVec.size(), 0.0);
//
//            for(int i=0; i<galSatVec.size(); ++i)
//            {
//                SatID sat( galSatVec[i] );
//
//                prefit(i) = gData.getValue(sat,TypeID::prefitC15);
//                weight(i) = gData.getValue(sat,TypeID::weight);
//            }
//
//            clock = dot(weight,prefit)/sum(weight);
//
//            for(int i=0; i<galSatVec.size(); ++i)
//            {
//                SatID sat( galSatVec[i] );
//
//                gData[sat][TypeID::cdtSta] = clock;
//            }

            return gData;
        }
        catch(Exception& u)
        {
            // Throw an exception if something unexpected happens
            ProcessingException e( getClassName() + u.what() );
            GPSTK_THROW(e);
        }

    }  // End of method 'ComputeStaClock::Process()'

}  // End of namespace gpstk
