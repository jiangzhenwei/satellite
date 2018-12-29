#pragma ident "$Id$"

/**
 * @file PhaseCodeAlignment.cpp
 * This class aligns phase with code measurements.
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
//  Dagoberto Salazar - gAGE ( http://www.gage.es ). 2008, 2009, 2011
//
//============================================================================


#include "PhaseCodeAlignment.hpp"

using namespace std;

namespace gpstk
{

    // Return a string identifying this object.
    std::string PhaseCodeAlignment::getClassName() const
    { return "PhaseCodeAlignment"; }


    /** Common constructor
     *
     * @param phase            Phase TypeID.
     * @param code             Code TypeID.
     * @param wavelength       Phase wavelength, in meters.
     * @param useArc           Whether satellite arcs will be used or not.
     */
    PhaseCodeAlignment::PhaseCodeAlignment( const SatID::SatelliteSystem& sys,
                                            const TypeID& phase,
                                            const TypeID& code,
                                            const double wavelength,
                                            bool useArc )
        : satSys(sys),
          phaseType(phase), codeType(code),
          useSatArcs(useArc), watchCSFlag(TypeID::CSL1)
    {
         // Set the wavelength
        setPhaseWavelength(wavelength);

    }  // End of 'PhaseCodeAlignment::PhaseCodeAlignment()'



    /** Method to set the phase wavelength to be used.
     *
     * @param wavelength       Phase wavelength, in meters.
     */
    PhaseCodeAlignment& PhaseCodeAlignment::setPhaseWavelength(double wavelength)
    {

         // Check that wavelength is bigger than zero
        if (wavelength > 0.0)
        {
            phaseWavelength = wavelength;
        }
        else
        {
            phaseWavelength = LC_WAVELENGTH_GPS_L1L2;   // Be default, GPS LC wavelength
        }

        return (*this);

    }  // End of 'PhaseCodeAlignment::setPhaseWavelength()'



    /** Return a satTypeValueMap object, adding the new data generated
     *  when calling this object.
     *
     * @param epoch     Time of observations.
     * @param gData     Data object holding the data.
     */
    satTypeValueMap& PhaseCodeAlignment::Process( const CommonTime& epoch,
                                                  satTypeValueMap& gData )
        throw(ProcessingException)
    {

        try
        {
            SatIDSet satRejectedSet;

            // Loop through all the satellites
            for(satTypeValueMap::iterator it = gData.begin();
                it != gData.end();
                ++it)
            {
                SatID sat( it->first );

                if(sat.system != satSys) continue;

                if( it->second.find(codeType)==it->second.end() ||
                    it->second.find(phaseType)==it->second.end() )
                    continue;

                // Check if satellite currently has entries
                SatData::const_iterator itDat( satData.find( sat ) );

                if( itDat == satData.end() )
                {
                    // If it doesn't have an entry, insert one
                    alignData aData;

                    satData[ sat ] = aData;
                }


                // Place to store if there was a cycle slip. False by default
                bool csflag(false);


                // Check if we want to use satellite arcs of cycle slip flags
                if(useSatArcs)
                {

                    double arcN(0.0);

                    try
                    {
                        // Try to extract the satellite arc value
                        arcN = it->second(TypeID::satArc);
                    }
                    catch(...)
                    {
                        // If satellite arc is missing, then schedule this
                        // satellite for removal
                        satRejectedSet.insert( sat );

                        continue;
                    }

                    // Check if satellite arc has changed
                    if( satData[sat].arcNumber != arcN )
                    {
                        // Set flag
                        csflag = true;

                        // Update satellite arc information
                        satData[sat].arcNumber = arcN;
                    }

                }  // End of first part of 'if(useSatArcs)'
                else
                {

                    double flag(0.0);

                    try
                    {
                        // Try to extract the CS flag value
                        flag = it->second(watchCSFlag);
                    }
                    catch(...)
                    {
                        // If flag is missing, then schedule this satellite
                        // for removal
                        satRejectedSet.insert( sat );

                        continue;
                    }

                    // Check if there was a cycle slip
                    if( flag > 0.0)
                    {
                        // Set flag
                        csflag = true;
                    }

                }  // End of second part of 'if(useSatArcs)...'


                // If there was an arc change or cycle slip, let's
                // compute the new offset
                if(csflag)
                {
                    // Compute difference between code and phase measurements
                    double diff( it->second(codeType) - it->second(phaseType) );

                    // Convert 'diff' to cycles
                    diff = diff/phaseWavelength;

                    // Convert 'diff' to an INTEGER number of cycles
                    diff = std::floor(diff);

                    // The new offset is the INTEGER number of cycles, in meters
                    satData[sat].offset = diff * phaseWavelength;
                }

                // Let's align the phase measurement using the
                // corresponding offset
                it->second[phaseType] = it->second[phaseType]
                                         + satData[sat].offset;


            }

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

    }  // End of 'PhaseCodeAlignment::Process()'



    /** Return a gnssSatTypeValue object, adding the new data generated
     *  when calling this object.
     *
     * @param gData    Data object holding the data.
     */
    gnssSatTypeValue& PhaseCodeAlignment::Process(gnssSatTypeValue& gData)
        throw(ProcessingException)
    {

        try
        {
            Process(gData.header.epoch, gData.body);

            return gData;
        }
        catch(Exception& u)
        {
            // Throw an exception if something unexpected happens
            ProcessingException e( getClassName() + ":" + u.what() );
            GPSTK_THROW(e);
        }

    }  // End of 'PhaseCodeAlignment::Process()'



    /** Return a gnssRinex object, adding the new data generated when
     *  calling this object.
     *
     * @param gData    Data object holding the data.
     */
    gnssRinex& PhaseCodeAlignment::Process(gnssRinex& gData)
        throw(ProcessingException)
    {

        try
        {
            Process(gData.header.epoch, gData.body);
            return gData;
        }
        catch(Exception& u)
        {
            // Throw an exception if something unexpected happens
            ProcessingException e( getClassName() + ":" + u.what() );
            GPSTK_THROW(e);
        }

    }  // End of 'PhaseCodeAlignment::Process()'


    /** Return a gnssDataMap object, adding the new data generated when
     *  calling this object.
     *
     * @param gData    Data object holding the data.
     */
    gnssDataMap& PhaseCodeAlignment::Process(gnssDataMap& gData)
        throw(ProcessingException)
    {

        for(gnssDataMap::iterator gdmIt = gData.begin();
            gdmIt != gData.end();
            ++gdmIt)
        {
            for(sourceDataMap::iterator sdmIt = gdmIt->second.begin();
                sdmIt != gdmIt->second.end();
                ++sdmIt)
            {
                SourceID source( sdmIt->first );

                SatDataMap::iterator satDataIt = satDataMap.find(source);

                if( satDataIt != satDataMap.end() )
                {
                    satData = satDataIt->second;
                }
                else
                {
                    satData = SatData();
                }

                Process( gdmIt->first, sdmIt->second );

                satDataMap[source] = satData;
            }
        }

        return gData;

    }  // End of 'PhaseCodeAlignment::Process()'


} // End of namespace gpstk
