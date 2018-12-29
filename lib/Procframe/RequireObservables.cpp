#pragma ident "$Id$"

/**
 * @file RequireObservables.cpp
 * This class filters out satellites with observations grossly out of bounds.
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


#include "RequireObservables.hpp"

namespace gpstk
{

    // Return a string identifying this object.
    std::string RequireObservables::getClassName() const
    { return "RequireObservables"; }



      /* Method to add a set of TypeID's to be required.
       *
       * @param typeSet    Set of TypeID's to be required.
       */
    RequireObservables& RequireObservables::addRequiredType( const SatID::SatelliteSystem& sys,
                                                             const TypeIDSet& typeSet )
    {
        if(sys == SatID::systemGPS)
        {
            requiredTypeSetOfGPS.insert( typeSet.begin(), typeSet.end() );
        }
        else if(sys == SatID::systemGLONASS)
        {
            requiredTypeSetOfGLO.insert( typeSet.begin(), typeSet.end() );
        }
        else if(sys == SatID::systemGalileo)
        {
            requiredTypeSetOfGAL.insert( typeSet.begin(), typeSet.end() );
        }
        else if(sys == SatID::systemBDS)
        {
            requiredTypeSetOfBDS.insert( typeSet.begin(), typeSet.end() );
        }

        return (*this);

    }  // End of method 'RequireObservables::addRequiredTypeOfGPS()'


      // Return a satTypeValueMap object, filtering the target observables.
      //
      // @param gData     Data object holding the data.
      //
    satTypeValueMap& RequireObservables::Process(satTypeValueMap& gData)
        throw(ProcessingException)
    {

        try
        {
            SatIDSet satRejectedSet;

            TypeIDSet requiredTypeSet;

            // Loop through all the satellites
            for(satTypeValueMap::iterator satIt = gData.begin();
                satIt != gData.end();
                ++satIt)
            {
                SatID sat( satIt->first );

                if(sat.system == SatID::systemGPS)
                {
                    requiredTypeSet = requiredTypeSetOfGPS;
                }
                else if(sat.system == SatID::systemGLONASS)
                {
                    requiredTypeSet = requiredTypeSetOfGLO;
                }
                else if(sat.system == SatID::systemGalileo)
                {
                    requiredTypeSet = requiredTypeSetOfGAL;
                }
                else if(sat.system == SatID::systemBDS)
                {
                    requiredTypeSet = requiredTypeSetOfBDS;
                }

                // Check all the indicated TypeID's
                if(requiredTypeSet.empty())
                {
                    satRejectedSet.insert( sat );
                    continue;
                }

                for(TypeIDSet::const_iterator typeIt = requiredTypeSet.begin();
                    typeIt != requiredTypeSet.end();
                    ++typeIt)
                {
                    // Try to find required type
                    typeValueMap::iterator it( satIt->second.find(*typeIt) );

                     // Now, check if this TypeID exists in this data structure
                    if(it == satIt->second.end())
                    {
                        // If we couldn't find type, then schedule this
                        // satellite for removal
                        satRejectedSet.insert( sat );

                        // It is not necessary to keep looking
                        //typeIt = requiredTypeSet.end();
                        //--typeIt;

                        break;
                    }
                }
            }

            // Let's remove satellites without all TypeID's
            gData.removeSatID(satRejectedSet);

            return gData;

        }
        catch(Exception& u)
        {
            // Throw an exception if something unexpected happens
            ProcessingException e( getClassName() + ":" + u.what() );
            GPSTK_THROW(e);
        }

    }  // End of 'RequireObservables::Process()'


      // Return a gnssDataMap object, filtering the target observables.
      //
      // @param gData     Data object holding the data.
      //
    gnssDataMap& RequireObservables::Process(gnssDataMap& gData)
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
                Process( sdmIt->second );
            }
        }

        return gData;

    }  // End of 'RequireObservables::Process()'


} // End of namespace gpstk
