#pragma ident "$Id$"

/**
 * @file CorrectObservables.cpp
 * This class corrects observables from effects such as antenna excentricity,
 * difference in phase centers, offsets due to tide effects, etc.
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
//  Dagoberto Salazar - gAGE ( http://www.gage.es ). 2007, 2008, 2009, 2011
//
//============================================================================


#include "CorrectObservables.hpp"

using namespace std;

namespace gpstk
{

    // Return a string identifying this object.
    std::string CorrectObservables::getClassName() const
    { return "CorrectObservables"; }



      /* Return a satTypeValueMap object, adding the new data generated
       * when calling this object.
       *
       * @param time      Epoch corresponding to the data.
       * @param gData     Data object holding the data.
       */
    satTypeValueMap& CorrectObservables::Process( const CommonTime& time,
                                                  satTypeValueMap& gData )
        throw(ProcessingException)
    {

        try
        {

            // Compute station latitude and longitude
            double lat( nominalPos.geodeticLatitude() );
            double lon( nominalPos.longitude() );

            // Define station position as a Triple, in ECEF
            Triple staPos( nominalPos.getX(),
                           nominalPos.getY(),
                           nominalPos.getZ() );

            // Compute initial displacement vectors, in meters [UEN]
            Triple initialBias( tide + monument );

            Triple dispL1( initialBias );
            Triple dispL2( initialBias );
            Triple dispL3( initialBias );
            Triple dispL5( initialBias );
            Triple dispL6( initialBias );
            Triple dispL7( initialBias );
            Triple dispL8( initialBias );
            Triple dispL9( initialBias );

            // Define a Triple that will hold satellite position, in ECEF
            Triple satPos(0.0, 0.0, 0.0);

            SatIDSet satRejectedSet;

            // Loop through all the satellites
            for(satTypeValueMap::iterator it = gData.begin();
                it != gData.end();
                ++it)
            {
                SatID sat( it->first );

                // Check if we have a valid Antenna object
                if( antenna.isValid() )
                {
                    // Compute phase center offsets
                    if( sat.system == SatID::systemGPS )
                    {
                        L1PhaseCenter = antenna.getAntennaEccentricity( Antenna::G01 );
                        L2PhaseCenter = antenna.getAntennaEccentricity( Antenna::G02 );
                        L5PhaseCenter = antenna.getAntennaEccentricity( Antenna::G02 );
                    }
                    else if( sat.system == SatID::systemGLONASS )
                    {
                        L1PhaseCenter = antenna.getAntennaEccentricity( Antenna::R01 );
                        L2PhaseCenter = antenna.getAntennaEccentricity( Antenna::R02 );
                    }
                    else if( sat.system == SatID::systemGalileo )
                    {
                        L1PhaseCenter = antenna.getAntennaEccentricity( Antenna::G01 );
                        L5PhaseCenter = antenna.getAntennaEccentricity( Antenna::G02 );
//                        L1PhaseCenter = antenna.getAntennaEccentricity( Antenna::E01 );
//                        L5PhaseCenter = antenna.getAntennaEccentricity( Antenna::E05 );
//                        L6PhaseCenter = antenna.getAntennaEccentricity( Antenna::E06 );
//                        L7PhaseCenter = antenna.getAntennaEccentricity( Antenna::E07 );
//                        L8PhaseCenter = antenna.getAntennaEccentricity( Antenna::E08 );
                    }
                    else if( sat.system == SatID::systemQZSS )
                    {
//                        L1PhaseCenter = antenna.getAntennaEccentricity( Antenna::J01 );
//                        L2PhaseCenter = antenna.getAntennaEccentricity( Antenna::J02 );
//                        L5PhaseCenter = antenna.getAntennaEccentricity( Antenna::J05 );
//                        L6PhaseCenter = antenna.getAntennaEccentricity( Antenna::J06 );
                    }
                    else if( sat.system == SatID::systemBDS )
                    {
//                        L2PhaseCenter = antenna.getAntennaEccentricity( Antenna::C02 );
//                        L6PhaseCenter = antenna.getAntennaEccentricity( Antenna::C06 );
//                        L7PhaseCenter = antenna.getAntennaEccentricity( Antenna::C07 );
                    }
                    else if( sat.system == SatID::systemIRNSS )
                    {
//                        L5PhaseCenter = antenna.getAntennaEccentricity( Antenna::I05 );
//                        L9PhaseCenter = antenna.getAntennaEccentricity( Antenna::I09 );
                    }
                }


                if( ( (*it).second.find(TypeID::satXECF) == (*it).second.end() ) ||
                    ( (*it).second.find(TypeID::satYECF) == (*it).second.end() ) ||
                    ( (*it).second.find(TypeID::satZECF) == (*it).second.end() ) )
                {
                    satRejectedSet.insert( sat );
                    continue;
                }
                else
                {
                    // Get satellite position out of GDS
                    satPos[0] = (*it).second[TypeID::satXECF];
                    satPos[1] = (*it).second[TypeID::satYECF];
                    satPos[2] = (*it).second[TypeID::satZECF];
                }


                // Declare the variables where antenna PC variations
                // will be stored.
                Triple L1Var( 0.0, 0.0, 0.0 );
                Triple L2Var( 0.0, 0.0, 0.0 );
                Triple L3Var( 0.0, 0.0, 0.0 );
                Triple L5Var( 0.0, 0.0, 0.0 );
                Triple L6Var( 0.0, 0.0, 0.0 );
                Triple L7Var( 0.0, 0.0, 0.0 );
                Triple L8Var( 0.0, 0.0, 0.0 );
                Triple L9Var( 0.0, 0.0, 0.0 );

                // Check if we have a valid Antenna object
                if( antenna.isValid() )
                {
                    // Check if we have elevation information
                    if( (*it).second.find(TypeID::elevation) != (*it).second.end() )
                    {

                        // Get elevation value
                        double elev( (*it).second[TypeID::elevation] );

                        // Check if azimuth is also required
                        if( !useAzimuth )
                        {
                            // In this case, use methods that only need elevation
                            try
                            {
                                // Compute phase center variation values
                                if( sat.system == SatID::systemGPS )
                                {
                                    L1Var = antenna.getAntennaPCVariation( Antenna::G01, elev );
                                    L2Var = antenna.getAntennaPCVariation( Antenna::G02, elev );
                                    L5Var = antenna.getAntennaPCVariation( Antenna::G02, elev );
                                }
                                else if( sat.system == SatID::systemGLONASS )
                                {
                                    L1Var = antenna.getAntennaPCVariation( Antenna::R01, elev );
                                    L2Var = antenna.getAntennaPCVariation( Antenna::R02, elev );
                                }
                                else if( sat.system == SatID::systemGalileo ) // no PCV available
                                {
                                    L1Var = antenna.getAntennaPCVariation( Antenna::G01, elev );
                                    L5Var = antenna.getAntennaPCVariation( Antenna::G02, elev );
//                                    L1Var = antenna.getAntennaPCVariation( Antenna::E01, elev );
//                                    L5Var = antenna.getAntennaPCVariation( Antenna::E05, elev );
//                                    L6Var = antenna.getAntennaPCVariation( Antenna::E06, elev );
//                                    L7Var = antenna.getAntennaPCVariation( Antenna::E07, elev );
//                                    L8Var = antenna.getAntennaPCVariation( Antenna::E08, elev );
                                }
                                else if( sat.system == SatID::systemQZSS ) // no PCV available
                                {
//                                    L1Var = antenna.getAntennaPCVariation( Antenna::J01, elev );
//                                    L2Var = antenna.getAntennaPCVariation( Antenna::J02, elev );
//                                    L5Var = antenna.getAntennaPCVariation( Antenna::J05, elev );
//                                    L6Var = antenna.getAntennaPCVariation( Antenna::J06, elev );
                                }
                                else if( sat.system == SatID::systemBDS ) // no PCV available
                                {
//                                    L2Var = antenna.getAntennaPCVariation( Antenna::C02, elev );
//                                    L6Var = antenna.getAntennaPCVariation( Antenna::C06, elev );
//                                    L7Var = antenna.getAntennaPCVariation( Antenna::C07, elev );
                                }
                                else if( sat.system == SatID::systemIRNSS ) // no PCV available
                                {
//                                    L5Var = antenna.getAntennaPCVariation( Antenna::I05, elev );
//                                    L9Var = antenna.getAntennaPCVariation( Antenna::I09, elev );
                                }
                            }
                            catch(InvalidRequest& ir)
                            {
                                // Throw an exception if something unexpected
                                // happens
                                ProcessingException e( getClassName() + ":"
                                        + "Unexpected problem found when trying to "
                                        + "compute antenna offsets" );

                                GPSTK_THROW(e);
                            }
                        }
                        else
                        {
                            // Check if we have azimuth information
                            if( (*it).second.find(TypeID::azimuth) !=
                                                            (*it).second.end() )
                            {

                                // Get azimuth value
                                double azim( (*it).second[TypeID::azimuth] );

                                // Use a gentle fallback mechanism to get antenna
                                // phase center variations
                                try
                                {
                                    // Compute phase center variation values
                                    if( sat.system == SatID::systemGPS )
                                    {
                                        L1Var = antenna.getAntennaPCVariation( Antenna::G01, elev, azim );
                                        L2Var = antenna.getAntennaPCVariation( Antenna::G02, elev, azim );
                                        L5Var = antenna.getAntennaPCVariation( Antenna::G02, elev, azim );
                                    }
                                    else if( sat.system == SatID::systemGLONASS )
                                    {
                                        L1Var = antenna.getAntennaPCVariation( Antenna::R01, elev, azim );
                                        L2Var = antenna.getAntennaPCVariation( Antenna::R02, elev, azim );
                                    }
                                    else if( sat.system == SatID::systemGalileo ) // no PCV available
                                    {
//                                        L1Var = antenna.getAntennaPCVariation( Antenna::E01, elev, azim );
//                                        L5Var = antenna.getAntennaPCVariation( Antenna::E05, elev, azim );
//                                        L6Var = antenna.getAntennaPCVariation( Antenna::E06, elev, azim );
//                                        L7Var = antenna.getAntennaPCVariation( Antenna::E07, elev, azim );
//                                        L8Var = antenna.getAntennaPCVariation( Antenna::E08, elev, azim );
                                    }
                                    else if( sat.system == SatID::systemQZSS ) // no PCV available
                                    {
//                                        L1Var = antenna.getAntennaPCVariation( Antenna::J01, elev, azim );
//                                        L2Var = antenna.getAntennaPCVariation( Antenna::J02, elev, azim );
//                                        L5Var = antenna.getAntennaPCVariation( Antenna::J05, elev, azim );
//                                        L6Var = antenna.getAntennaPCVariation( Antenna::J06, elev, azim );
                                    }
                                    else if( sat.system == SatID::systemBDS ) // no PCV available
                                    {
//                                        L2Var = antenna.getAntennaPCVariation( Antenna::C02, elev, azim );
//                                        L6Var = antenna.getAntennaPCVariation( Antenna::C06, elev, azim );
//                                        L7Var = antenna.getAntennaPCVariation( Antenna::C07, elev, azim );
                                    }
                                    else if( sat.system == SatID::systemIRNSS ) // no PCV available
                                    {
//                                        L5Var = antenna.getAntennaPCVariation( Antenna::I05, elev, azim );
//                                        L9Var = antenna.getAntennaPCVariation( Antenna::I09, elev, azim );
                                    }
                                }
                                catch(InvalidRequest& ir)
                                {
                                    // We  "graceful degrade" to a simpler mechanism
                                    try
                                    {
                                        // Compute phase center variation values
                                        if( sat.system == SatID::systemGPS )
                                        {
                                            L1Var = antenna.getAntennaPCVariation( Antenna::G01, elev );
                                            L2Var = antenna.getAntennaPCVariation( Antenna::G02, elev );
                                            L5Var = antenna.getAntennaPCVariation( Antenna::G02, elev );
                                        }
                                        else if( sat.system == SatID::systemGLONASS )
                                        {
                                            L1Var = antenna.getAntennaPCVariation( Antenna::R01, elev );
                                            L2Var = antenna.getAntennaPCVariation( Antenna::R02, elev );
                                        }
                                        else if( sat.system == SatID::systemGalileo )
                                        {
//                                            L1Var = antenna.getAntennaPCVariation( Antenna::E01, elev );
//                                            L5Var = antenna.getAntennaPCVariation( Antenna::E05, elev );
//                                            L6Var = antenna.getAntennaPCVariation( Antenna::E06, elev );
//                                            L7Var = antenna.getAntennaPCVariation( Antenna::E07, elev );
//                                            L8Var = antenna.getAntennaPCVariation( Antenna::E08, elev );
                                        }
                                        else if( sat.system == SatID::systemQZSS )
                                        {
//                                            L1Var = antenna.getAntennaPCVariation( Antenna::J01, elev );
//                                            L2Var = antenna.getAntennaPCVariation( Antenna::J02, elev );
//                                            L5Var = antenna.getAntennaPCVariation( Antenna::J05, elev );
//                                            L6Var = antenna.getAntennaPCVariation( Antenna::J06, elev );
                                        }
                                        else if( sat.system == SatID::systemBDS )
                                        {
//                                            L1Var = antenna.getAntennaPCVariation( Antenna::C01, elev );
//                                            L2Var = antenna.getAntennaPCVariation( Antenna::C02, elev );
                                        }
                                        else if( sat.system == SatID::systemIRNSS )
                                        {
//                                            L5Var = antenna.getAntennaPCVariation( Antenna::I05, elev );
//                                            L9Var = antenna.getAntennaPCVariation( Antenna::I09, elev );
                                        }
                                    }
                                    catch(InvalidRequest& ir)
                                    {
                                        // Throw an exception if something unexpected
                                        // happens
                                        ProcessingException e( getClassName() + ":"
                                            + "Unexpected problem found when trying to "
                                            + "compute antenna offsets" );

                                        GPSTK_THROW(e);
                                    }
                                }
                            }
                            else
                            {
                                // Throw an exception if something unexpected happens
                                ProcessingException e( getClassName() + ":"
                                    + "Azimuth information could not be found, "
                                    + "so antenna PC offsets can not be computed");

                                GPSTK_THROW(e);

                            }  // End of 'if( (*it).second.find(TypeID::azimuth) !=...'

                        }  // End of 'if( !useAzimuth )'
                    }
                    else
                    {
                        // Throw an exception if there is no elevation data
                        ProcessingException e( getClassName() + ":"
                            + "Elevation information could not be found, "
                            + "so antenna PC offsets can not be computed" );

                        GPSTK_THROW(e);

                    }  // End of 'if( (*it).second.find(TypeID::elevation) != ...'

                }  // End of 'if( antenna.isValid() )...'


                // Update displacement vectors with current phase centers
                Triple dL1( dispL1 + L1PhaseCenter );
                Triple dL2( dispL2 + L2PhaseCenter );
                Triple dL3( dispL3 + L3PhaseCenter );
                Triple dL5( dispL5 + L5PhaseCenter );
                Triple dL6( dispL6 + L6PhaseCenter );
                Triple dL7( dispL7 + L7PhaseCenter );
                Triple dL8( dispL8 + L8PhaseCenter );
                Triple dL9( dispL9 + L9PhaseCenter );

                // Compute vector station-satellite, in ECF
                Triple ray(staPos - satPos);

                ray = ray.unitVector();

                // Rotate vector ray to UEN reference frame
                ray = (ray.R3(lon)).R2(-lat);

                // Convert ray to an unitary vector
                ray = ray.unitVector();

                // Compute corrections = displacement vectors components
                // along ray direction.
                double corrL1( dL1.dot(ray) + L1Var[0] );
                double corrL2( dL2.dot(ray) + L2Var[0] );
                double corrL3( dL3.dot(ray) + L3Var[0] );
                double corrL5( dL5.dot(ray) + L5Var[0] );
                double corrL6( dL6.dot(ray) + L6Var[0] );
                double corrL7( dL7.dot(ray) + L7Var[0] );
                double corrL8( dL8.dot(ray) + L8Var[0] );
                double corrL9( dL9.dot(ray) + L9Var[0] );


                // Find which observables are present, and then
                // apply corrections

                // Look for C1
                if( (*it).second.find(TypeID::C1) != (*it).second.end() )
                {
                   (*it).second[TypeID::C1] = (*it).second[TypeID::C1] - corrL1;
                }

                // Look for L1
                if( (*it).second.find(TypeID::L1) != (*it).second.end() )
                {
                   (*it).second[TypeID::L1] = (*it).second[TypeID::L1] - corrL1;
                }

                // Look for C2
                if( (*it).second.find(TypeID::C2) != (*it).second.end() )
                {
                   (*it).second[TypeID::C2] = (*it).second[TypeID::C2] - corrL2;
                }

                // Look for L2
                if( (*it).second.find(TypeID::L2) != (*it).second.end() )
                {
                   (*it).second[TypeID::L2] = (*it).second[TypeID::L2] - corrL2;
                }

                // Look for C3
                if( (*it).second.find(TypeID::C3) != (*it).second.end() )
                {
                   (*it).second[TypeID::C3] = (*it).second[TypeID::C3] - corrL3;
                }

                // Look for L3
                if( (*it).second.find(TypeID::L3) != (*it).second.end() )
                {
                   (*it).second[TypeID::L3] = (*it).second[TypeID::L3] - corrL3;
                }

                // Look for C5
                if( (*it).second.find(TypeID::C5) != (*it).second.end() )
                {
                   (*it).second[TypeID::C5] = (*it).second[TypeID::C5] - corrL5;
                }

                // Look for L5
                if( (*it).second.find(TypeID::L5) != (*it).second.end() )
                {
                   (*it).second[TypeID::L5] = (*it).second[TypeID::L5] - corrL5;
                }

                // Look for C6
                if( (*it).second.find(TypeID::C6) != (*it).second.end() )
                {
                   (*it).second[TypeID::C6] = (*it).second[TypeID::C6] - corrL6;
                }

                // Look for L6
                if( (*it).second.find(TypeID::L6) != (*it).second.end() )
                {
                   (*it).second[TypeID::L6] = (*it).second[TypeID::L6] - corrL6;
                }

                // Look for C7
                if( (*it).second.find(TypeID::C7) != (*it).second.end() )
                {
                   (*it).second[TypeID::C7] = (*it).second[TypeID::C7] - corrL7;
                }

                // Look for L7
                if( (*it).second.find(TypeID::L7) != (*it).second.end() )
                {
                   (*it).second[TypeID::L7] = (*it).second[TypeID::L7] - corrL7;
                }

                // Look for C8
                if( (*it).second.find(TypeID::C8) != (*it).second.end() )
                {
                   (*it).second[TypeID::C8] = (*it).second[TypeID::C8] - corrL8;
                }

                // Look for L8
                if( (*it).second.find(TypeID::L8) != (*it).second.end() )
                {
                   (*it).second[TypeID::L8] = (*it).second[TypeID::L8] - corrL8;
                }

                // Look for C9
                if( (*it).second.find(TypeID::C9) != (*it).second.end() )
                {
                   (*it).second[TypeID::C9] = (*it).second[TypeID::C9] - corrL9;
                }

                // Look for L9
                if( (*it).second.find(TypeID::L9) != (*it).second.end() )
                {
                   (*it).second[TypeID::L9] = (*it).second[TypeID::L9] - corrL9;
                }

            } // End of 'for(satTypeValueMap::iterator = it, ...)'


            // Remove satellites with missing data
            gData.removeSatID(satRejectedSet);

            return gData;

        }
        catch(Exception& u)
        {
            // Throw an exception if something unexpected happens
            ProcessingException e( getClassName() + ":" + u.what() );
            GPSTK_THROW(e);
        }

    }  // End of method 'CorrectObservables::Process()'


     /** Return a gnssDataMap object, adding the new data generated when
      *  calling this object.
      *
      * @param gData    Data object holding the data.
      */
    gnssDataMap& CorrectObservables::Process(gnssDataMap& gData)
        throw(ProcessingException)
    {

        SourceIDSet sourceRejectedSet;

        for( gnssDataMap::iterator gdmIt = gData.begin();
             gdmIt != gData.end();
             ++gdmIt )
        {
            CommonTime epoch( gdmIt->first );
            epoch.setTimeSystem( TimeSystem::Unknown );

            for( sourceDataMap::iterator sdmIt = gdmIt->second.begin();
                 sdmIt != gdmIt->second.end();
                 ++sdmIt )
            {
                SourceID source( sdmIt->first );
                string station( source.sourceName.substr(0,4) );

                string upper( station );
                transform(upper.begin(),upper.end(),upper.begin(),::toupper);

                SourceID temp( SourceID(SourceID::Mixed, upper) );

                if( sourceMonument.find(temp) == sourceMonument.end() ||
                    sourceAntenna.find(temp) == sourceAntenna.end() )
                {
                    sourceRejectedSet.insert( source );
                    continue;
                }

                monument = sourceMonument[temp];
                antenna = sourceAntenna[temp];

                if(pMSCStore == NULL)
                {
                    sourceRejectedSet.insert( source );
                    continue;
                }

                MSCData mscData;
                try
                {
                    mscData = pMSCStore->findMSC(upper,epoch);
                }
                catch(...)
                {
                    sourceRejectedSet.insert( source );
                    continue;
                }

                nominalPos = mscData.coordinates;

                pStaTides->setNominalPosition( nominalPos );
                tide = pStaTides->getTides( epoch, upper );

                Process( gdmIt->first, sdmIt->second );
            }
        }

        gData.removeSourceID( sourceRejectedSet );

        return gData;

    } // End of method 'CorrectObservables::Process()'


}  // End of namespace gpstk
