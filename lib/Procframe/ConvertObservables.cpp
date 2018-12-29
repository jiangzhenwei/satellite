#pragma ident "$Id$"

/**
* @file ConvertObservables.cpp
* Class to convert observables, from 3-char to 2-char.
*/

#include "ConvertObservables.hpp"

using namespace std;

namespace gpstk
{

      // Return a string identifying this object.
    string ConvertObservables::getClassName() const
    { return "ConvertObservables"; }


      /* Return a satTypeValueMap object, adding the new data generated
       *  when calling this object.
       *
       * @param time      Epoch corresponding to the data.
       * @param gData     Data object holding the data.
       */
    satTypeValueMap& ConvertObservables::Process( const CommonTime& time,
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

                typeValueMap::iterator itType;

                if(sat.system == SatID::systemGPS)
                {
                    // C1
                    itType = it->second.find(TypeID::C1W);
                    if( itType != it->second.end() )
                    {
                        it->second[TypeID::C1] = it->second[TypeID::C1W];
                    }
                    else
                    {
//                        satRejectedSet.insert( sat );
                    }

                    // L1
                    itType = it->second.find(TypeID::L1C);
                    if( itType != it->second.end() )
                    {
                        it->second[TypeID::L1] = it->second[TypeID::L1C];
                    }
                    else
                    {
//                        satRejectedSet.insert( sat );
                    }

                    // C2
                    itType = it->second.find(TypeID::C2W);
                    if( itType != it->second.end() )
                    {
                        it->second[TypeID::C2] = it->second[TypeID::C2W];
                    }
                    else
                    {
//                        satRejectedSet.insert( sat );
                    }

                    // L2
                    itType = it->second.find(TypeID::L2W);
                    if( itType != it->second.end() )
                    {
                        it->second[TypeID::L2] = it->second[TypeID::L2W];
                    }
                    else
                    {
//                        satRejectedSet.insert( sat );
                    }

                    // C5
                    itType = it->second.find(TypeID::C5X);
                    if( itType != it->second.end() )
                    {
                        it->second[TypeID::C5] = it->second[TypeID::C5X];
                    }
                    else
                    {
//                        satRejectedSet.insert( sat );
                    }

                    // L5
                    itType = it->second.find(TypeID::L5X);
                    if( itType != it->second.end() )
                    {
                        it->second[TypeID::L5] = it->second[TypeID::L5X];
                    }
                    else
                    {
//                        satRejectedSet.insert( sat );
                    }
                }
                else if(sat.system == SatID::systemGLONASS)
                {
                    satRejectedSet.insert(sat);
                }
                else if(sat.system == SatID::systemGalileo)
                {
                    // C1
                    itType = it->second.find(TypeID::C1X);
                    if( itType != it->second.end() )
                    {
                        it->second[TypeID::C1] = it->second[TypeID::C1X];
                    }
                    else
                    {
//                        satRejectedSet.insert( sat );
                    }

                    // L1
                    itType = it->second.find(TypeID::L1X);
                    if( itType != it->second.end() )
                    {
                        it->second[TypeID::L1] = it->second[TypeID::L1X];
                    }
                    else
                    {
                        satRejectedSet.insert( sat );
                    }

                    // C5
                    itType = it->second.find(TypeID::C5X);
                    if( itType != it->second.end() )
                    {
                        it->second[TypeID::C5] = it->second[TypeID::C5X];
                    }
                    else
                    {
//                        satRejectedSet.insert( sat );
                    }

                    // L5
                    itType = it->second.find(TypeID::L5X);
                    if( itType != it->second.end() )
                    {
                        it->second[TypeID::L5] = it->second[TypeID::L5X];
                    }
                    else
                    {
//                        satRejectedSet.insert( sat );
                    }
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

    }  // End of method 'ConvertObservables::Process()'


     /** Return a gnssDataMap object, adding the new data generated when
      *  calling a modeling object.
      *
      * @param gData    Data object holding the data.
      */
    gnssDataMap& ConvertObservables::Process(gnssDataMap& gData)
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

    }  // End of method 'ConvertObservables::Process()'

}  // End of namespace gpstk
