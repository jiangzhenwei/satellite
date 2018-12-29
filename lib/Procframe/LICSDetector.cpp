#pragma ident "$Id$"

/**
 * @file LICSDetector.cpp
 * This is a class to detect cycle slips using LI observables.
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
//  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110, USA
//
//  Dagoberto Salazar - gAGE ( http://www.gage.es ). 2007, 2008, 2011
//
//============================================================================


#include "LICSDetector.hpp"


namespace gpstk
{

      // Return a string identifying this object.
    std::string LICSDetector::getClassName() const
    { return "LICSDetector"; }


      /* Common constructor
       *
       * @param mThr    Minimum threshold to declare cycle slip, in meters.
       * @param drift   LI combination limit drift, in meters/second.
       * @param dtMax   Maximum interval of time allowed between two
       *                successive epochs, in seconds.
       */
    LICSDetector::LICSDetector( const double& mThr,
                                const double& drift,
                                const double& dtMax,
                                const bool& use )
        : obsType(TypeID::LI12),
          lliType1(TypeID::LLI1), lliType2(TypeID::LLI2),
          resultType1(TypeID::CSL1), resultType2(TypeID::CSL2),
          useLLI(use)
    {
        setDeltaTMax(dtMax);
        setMinThreshold(mThr);
        setLIDrift(drift);
    }


      /* Return a satTypeValueMap object, adding the new data generated
       *  when calling this object.
       *
       * @param epoch     Time of observations.
       * @param gData     Data object holding the data.
       * @param epochflag Epoch flag.
       */
    satTypeValueMap& LICSDetector::Process( const CommonTime& epoch,
                                            satTypeValueMap& gData,
                                            const short& epochflag )
        throw(ProcessingException)
    {

        try
        {
            double value1(0.0);
            double lli1(0.0);
            double lli2(0.0);

            SatIDSet satRejectedSet;

            // Loop through all the satellites
            for(satTypeValueMap::iterator it = gData.begin();
                it != gData.end();
                ++it)
            {
                SatID sat = (*it).first;

                if(sat.system != satSys) continue;

                try
                {
                    // Try to extract the values
                    value1 = (*it).second(obsType);
                }
                catch(...)
                {
                    // If some value is missing, then schedule this satellite
                    // for removal
                    satRejectedSet.insert( sat );
                    continue;
                }

                if (useLLI)
                {
                    try
                    {
                        // Try to get the LLI1 index
                        lli1  = (*it).second(lliType1);
                    }
                    catch(...)
                    {
                        // If LLI #1 is not found, set it to zero
                        // You REALLY want to have BOTH LLI indexes properly set
                        lli1 = 0.0;
                    }

                    try
                    {
                        // Try to get the LLI2 index
                        lli2  = (*it).second(lliType2);
                    }
                    catch(...)
                    {
                        // If LLI #2 is not found, set it to zero
                        // You REALLY want to have BOTH LLI indexes properly set
                        lli2 = 0.0;
                    }
                }

                // If everything is OK, then get the new values inside the
                // structure. This way of computing it allows concatenation of
                // several different cycle slip detectors
                (*it).second[resultType1] += getDetection( epoch,
                                                           sat,
                                                           (*it).second,
                                                           epochflag,
                                                           value1,
                                                           lli1,
                                                           lli2 );

                if ( (*it).second[resultType1] > 1.0 )
                {
                    (*it).second[resultType1] = 1.0;
                }

                // We will mark both cycle slip flags
                (*it).second[resultType2] = (*it).second[resultType1];

            }

            // Remove satellites with missing data
            if(checkExist)
            {
                gData.removeSatID(satRejectedSet);
            }

            return gData;
        }
        catch(Exception& u)
        {
            // Throw an exception if something unexpected happens
            ProcessingException e( getClassName() + ":" + u.what() );
            GPSTK_THROW(e);
        }

    }  // End of method 'LICSDetector::Process()'


      /* Return a gnssRinex object, adding the new data generated when
       * calling this object.
       *
       * @param gData    Data object holding the data.
       */
    gnssRinex& LICSDetector::Process(gnssRinex& gData)
        throw(ProcessingException)
    {

        try
        {
            Process(gData.header.epoch, gData.body, gData.header.epochFlag);

            return gData;
        }
        catch(Exception& u)
        {
            // Throw an exception if something unexpected happens
            ProcessingException e( getClassName() + ":" + u.what() );
            GPSTK_THROW(e);
        }

    }  // End of method 'LICSDetector::Process()'


      /* Return a gnssDataMap object, adding the new data generated when
       * calling this object.
       *
       * @param gData    Data object holding the data.
       */
    gnssDataMap& LICSDetector::Process(gnssDataMap& gData)
        throw(ProcessingException)
    {

        SourceID source;

        for(gnssDataMap::iterator gdmIt = gData.begin();
            gdmIt != gData.end();
            ++gdmIt)
        {
            for(sourceDataMap::iterator sdmIt = gdmIt->second.begin();
                sdmIt != gdmIt->second.end();
                ++sdmIt)
            {
                source = sdmIt->first;

                LIDataMap::iterator liDataMapIt = liDataMap.find(source);

                if( liDataMapIt != liDataMap.end() )
                {
                    liData = liDataMapIt->second;
                }
                else
                {
                    liData = LIData();
                }

                Process( gdmIt->first, sdmIt->second );

                liDataMap[source] = liData;
            }
       }

       return gData;

    }  // End of method 'LICSDetector::Process()'


      /* Method that implements the LI cycle slip detection algorithm
       *
       * @param epoch     Time of observations.
       * @param sat       SatID.
       * @param tvMap     Data structure of TypeID and values.
       * @param epochflag Epoch flag.
       * @param li        Current LI observation value.
       * @param lli1      LLI1 index.
       * @param lli2      LLI2 index.
       */
    double LICSDetector::getDetection( const CommonTime& epoch,
                                       const SatID& sat,
                                       typeValueMap& tvMap,
                                       const short& epochflag,
                                       const double& li,
                                       const double& lli1,
                                       const double& lli2 )
    {

        bool reportCS(false);

         // Difference between current and former epochs, in sec
        double currentDeltaT(0.0);

         // Difference between current and former LI values
        double currentBias(0.0);

         // Limit to declare cycle slip
        double deltaLimit(0.0);

        double delta(0.0);
        double tempLLI1(0.0);
        double tempLLI2(0.0);


         // Get the difference between current epoch and former epoch,
         // in seconds
        currentDeltaT = ( epoch - liData[sat].formerEpoch );

         // Store current epoch as former epoch
        liData[sat].formerEpoch = epoch;

         // Current value of LI difference
        currentBias = li - liData[sat].formerLI;

         // Increment window size
        ++liData[sat].windowSize;

         // Check if receiver already declared cycle slip or too much time
         // has elapsed
         // Note: If tvMap(lliType1) or tvMap(lliType2) don't exist, then 0
         // will be returned and those tests will pass
        if( (tvMap(lliType1) == 1.0) ||
            (tvMap(lliType1) == 3.0) ||
            (tvMap(lliType1) == 5.0) ||
            (tvMap(lliType1) == 7.0) )
        {
            tempLLI1 = 1.0;
        }

        if( (tvMap(lliType2) == 1.0) ||
            (tvMap(lliType2) == 3.0) ||
            (tvMap(lliType2) == 5.0) ||
            (tvMap(lliType2) == 7.0) )
        {
            tempLLI2 = 1.0;
        }

//        if( (epochflag==1)  ||
//            (epochflag==6)  ||
//            (tempLLI1==1.0) ||
//            (tempLLI2==1.0) ||
//            (currentDeltaT > deltaTMax) )
         /**
          *  The 'epochflag' is not reliable.
          */
        if( (tempLLI1==1.0) ||
            (tempLLI2==1.0) ||
            (currentDeltaT > deltaTMax) )
        {
            // We reset the filter with this
            liData[sat].windowSize = 0;
            reportCS = true;
        }

        if(liData[sat].windowSize > 1)
        {
            deltaLimit = minThreshold + std::abs(LIDrift*currentDeltaT);

            // Compute a linear interpolation and compute
            // LI_predicted - LI_current
            delta = std::abs( currentBias - (liData[sat].formerBias *
                            currentDeltaT / liData[sat].formerDeltaT) );

            if (delta > deltaLimit)
            {
                // We reset the filter with this
                liData[sat].windowSize = 0;
                reportCS = true;
            }
        }

         // Let's prepare for the next time
        liData[sat].formerLI = li;
        liData[sat].formerBias = currentBias;
        liData[sat].formerDeltaT = currentDeltaT;

        if(reportCS)
        {
            return 1.0;
        }
        else
        {
            return 0.0;
        }

    }  // End of method 'LICSDetector::getDetection()'


}  // End of namespace gpstk
