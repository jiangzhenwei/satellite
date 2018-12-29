#pragma ident "$Id$"

/**
 * @file ComputeLinear.cpp
 * This class computes linear combinations of GDS data.
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


#include "ComputeLinear.hpp"

using namespace std;

namespace gpstk
{

      // Return a string identifying this object.
    std::string ComputeLinear::getClassName() const
    { return "ComputeLinear"; }



      /* Return a satTypeValueMap object, adding the new data generated when
       * calling this object.
       *
       * @param time      Epoch corresponding to the data.
       * @param gData     Data object holding the data.
       */
    satTypeValueMap& ComputeLinear::Process( const CommonTime& time,
                                             satTypeValueMap& gData )
        throw(ProcessingException)
    {

        try
        {
            LinearCombList linearList;

            // Loop through all the satellites
            for(satTypeValueMap::iterator it = gData.begin();
                it != gData.end();
                ++it)
            {
                SatID sat = (*it).first;

                if(sat.system == SatID::systemGPS)
                    linearList = linearListOfGPS;
                else if(sat.system == SatID::systemGalileo)
                    linearList = linearListOfGAL;
                else if(sat.system == SatID::systemBDS)
                    linearList = linearListOfBDS;

                // Loop through all the defined linear combinations
                for(LinearCombList::const_iterator pos = linearList.begin();
                    pos != linearList.end();
                    ++pos)
                {
                    double result(0.0);

                    bool exist( true );

                    // Read the information of each linear combination
                    for(typeValueMap::const_iterator iter = pos->body.begin();
                        iter != pos->body.end();
                        ++iter)
                    {
                        double temp(0.0);

                        TypeID type(iter->first);

                        if( (*it).second.find(type) != (*it).second.end() )
                        {
                            temp = (*it).second[type];
                        }
                        else
                        {
                            exist = false;
                            temp = 0.0;
                        }

                        result = result + (*iter).second * temp;
                    }

                    // Store the result in the proper place
                    if(!checkExist || exist)
                    {
                        (*it).second[pos->header] = result;
                    }

                }

            }

            return gData;

        }
        catch(Exception& u)
        {
            // Throw an exception if something unexpected happens
            ProcessingException e( getClassName() + ":" + u.what() );
            GPSTK_THROW(e);
        }

    }  // End of method 'ComputeLinear::Process()'



     /** Return a gnssDataMap object, adding the new data generated
      *  when calling this object.
      *
      * @param gData    Data object holding the data.
      */
    gnssDataMap& ComputeLinear::Process(gnssDataMap& gData)
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
                Process( gdmIt->first, sdmIt->second );
            }
        }

        return gData;

    }  // End of method 'ComputeLinear::Process()'


} // End of namespace gpstk
