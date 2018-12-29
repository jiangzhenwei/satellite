#pragma ident "$Id$"

/**
* @file GPSTK_CC2NONCC_CPP
* Class to convert CC(cross-correlation) to NONCC(non cross-correlation).
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
//  Copyright (c)
//
//  Q.Liu, Wuhan Uniersity, 2015
//============================================================================
//  Modification
//  2016-1-26    For C1/P2 receiver, the observable list maybe has P1(0),then
//               this Satelite would be filtered, this result is unreasonable,
//               so we delete the P1 filter.(Q.Liu)
//
//  2016-4-18    Add the exception when DCB file is not provided, and we will
//               not do the DCB correction and will output the error message,
//               but not stop the data processing work.(Q.Liu)
//============================================================================

#include "CC2NONCC.hpp"


using namespace std;

namespace gpstk
{

      // Return a string identifying this object.
    string CC2NONCC::getClassName() const
    { return "CC2NONCC"; }


      /* Return a satTypeValueMap object, adding the new data generated
       *  when calling this object.
       *
       * @param time      Epoch corresponding to the data.
       * @param gData     Data object holding the data.
       */
    satTypeValueMap& CC2NONCC::Process( const CommonTime& time,
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

                if(sat.system == SatID::systemGPS)
                {
                    typeValueMap::iterator itC1W( it->second.find(TypeID::C1W) );

                    if(itC1W != it->second.end())
                    {
                        continue;
                    }

                    if(dcbP1C1 != NULL)
                    {
                        // Get the Sat's DCB value
                        double Bp1c1(0.0);      // in ns
                        try
                        {
                            Bp1c1 = (*dcbP1C1).getDCB(sat);
                        }
                        catch(...)
                        {
                            Bp1c1 = 0.0;
                        }

                        typeValueMap::iterator itC1C = it->second.find(TypeID::C1C);

                        if(itC1C == it->second.end())
                        {
                            satRejectedSet.insert( sat );
                            continue;
                        }

                        // For receiver noncc (C1,P2)
                        // For the noncc: only C1 should be corrected
                        // C1 -> C1 + (P1-C1)
                        if(itC1C != it->second.end())
                        {
                            // Correct
                            it->second[TypeID::C1W] = it->second[TypeID::C1C] +
                                                      Bp1c1 * C_MPS * 1.0e-9;
                        }
                        continue;

                    } // End of 'if(dcbP1C1 != NULL)'
                }
                else if(sat.system == SatID::systemGLONASS)
                {
                    satRejectedSet.insert(sat);
                    continue;
                }
                else if(sat.system == SatID::systemGalileo)
                {
                    continue;
                }
                else if(sat.system == SatID::systemBDS)
                {
                    continue;
                }
                else
                {
                    satRejectedSet.insert(sat);
                    continue;
                }

            }  // End of 'for (it = gData.begin(); it != gData.end(); ++it)'

            gData.removeSatID(satRejectedSet);

            return gData;
        }
        catch(Exception& u)
        {
            // Throw an exception if something unexpected happens
            ProcessingException e( getClassName() + ":" + u.what() );
            GPSTK_THROW(e);
        }

    }  // End of method 'CC2NONCC::Process()'


    /** Return a gnssDataMap object, adding the new data generated
     *  when calling this object.
     *
     * @param gData     Data object holding the data.
     */
    gnssDataMap& CC2NONCC::Process(gnssDataMap& gData)
        throw(ProcessingException)
    {
        SourceIDSet sourceRejectedSet;

        for( gnssDataMap::iterator gdmIt = gData.begin();
             gdmIt != gData.end();
             ++gdmIt )
        {
            for( sourceDataMap::iterator sdmIt = gdmIt->second.begin();
                 sdmIt != gdmIt->second.end();
                 ++sdmIt )
            {
                Process( gdmIt->first, sdmIt->second );
            }
        }

        gData.removeSourceID( sourceRejectedSet );

        return gData;

    }  // End of method 'CC2NONCC::Process()'

}  // End of namespace gpstk
