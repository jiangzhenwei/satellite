#pragma ident "$Id$"

/**
 * @file MWFilter.cpp
 * This is a class to filter MW observations.
 */


#include "MWFilter.hpp"


using namespace std;


namespace gpstk
{


      // Return a string identifying this object.
    std::string MWFilter::getClassName() const
    { return "MWFilter"; }



      /* Return a satTypeValueMap object, adding the new data generated
       * when calling this object.
       *
       * @param epoch     Time of observations.
       * @param gData     Data object holding the data.
       */
    satTypeValueMap& MWFilter::Process( const CommonTime& epoch,
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
                SatID sat( (*it).first );

                double mw(0.0);

                try
                {
                    // Try to extract the values
                    mw = (*it).second(TypeID::MW12);
                }
                catch(...)
                {
                    // If some value is missing, then schedule this satellite
                    // for removal
                    satRejectedSet.insert( sat );
                    continue;
                }

                // First check if this satellite has previous arc information
                if( arcData.find( sat ) == arcData.end() )
                {
                    // If it doesn't have an entry, insert one
                    arcData[ sat ] = 0.0;
                }

                //>>> Start to filter the MW data

                // Then, check both if there is arc information, and if current
                // arc number is different from arc number in storage (which
                // means a cycle slip happened), reset the buffer for smoothing data.
                if( (*it).second.find(TypeID::satArc) != (*it).second.end() &&
                    (*it).second(TypeID::satArc) != arcData[ sat ] )
                {
                    // If different, update satellite arc in storage
                    arcData[ sat ] = (*it).second(TypeID::satArc);

                    // We reset the filter with this
                    mwData[ sat ].windowSize = 1;

                    // Now, the mean value equals with the current 'ambValue'
                    mwData[ sat ].meanMW = mw;

                    // Now, the variance value equals with the current 'ambValue'
                    mwData[ sat ].varMW = 0.25*0.25;
                }
                else
                {
                    // Increment window size
                    ++mwData[sat].windowSize;

                    // MW bias from the mean value
                    double mwBias(mw - mwData[sat].meanMW);
                    double size( static_cast<double>(mwData[sat].windowSize) );

                    // Compute average
                    // meanMW(i)= meanMW(i-1) + ( mwBias ) / size;
                    mwData[sat].meanMW += mwBias / size;

                    // Compute variance
                    // Var(i) = Var(i-1) + [ ( mw(i) - meanMW)^2 - Var(i-1) ]/(i);
                    mwData[sat].varMW  += ( mwBias*mwBias - mwData[sat].varMW ) / size;

                }  // End of smoothing data

                // Insert the mean value into gData
                (*it).second[TypeID::meanMW12] = mwData[sat].meanMW;
                (*it).second[TypeID::varMW12] = mwData[sat].varMW;

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

    }  // End of method 'MWFilter::Process()'



      /* Return a gnssRinex object, adding the new data generated when
       * calling this object.
       *
       * @param gData    Data object holding the data.
       */
    gnssRinex& MWFilter::Process(gnssRinex& gData)
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

    }  // End of method 'MWFilter::Process()'


      /* Return a gnssDataMap object, adding the new data generated when
       * calling this object.
       *
       * @param gData    Data object holding the data.
       */
    gnssDataMap& MWFilter::Process(gnssDataMap& gData)
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

                // arc data
                ArcDataMap::iterator arcDataIt = arcDataMap.find(source);

                if(arcDataIt != arcDataMap.end())
                {
                    arcData = arcDataIt->second;
                }
                else
                {
                    arcData = ArcData();
                }

                // mw data
                MWDataMap::iterator mwDataIt = mwDataMap.find(source);

                if(mwDataIt != mwDataMap.end())
                {
                    mwData = mwDataIt->second;
                }
                else
                {
                    mwData = MWData();
                }

                Process( gdmIt->first, sdmIt->second );

                arcDataMap[source] = arcData;
                mwDataMap[source] = mwData;
            }
        }

        return gData;

    }  // End of method 'MWFilter::Process()'


}  // End of namespace gpstk
