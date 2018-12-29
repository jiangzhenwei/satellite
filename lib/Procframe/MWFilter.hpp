#pragma ident "$Id$"

/**
 * @file MWFilter.hpp
 * This is a class to filter MW observations.
 */

#ifndef GPSTK_MWFILTER_HPP
#define GPSTK_MWFILTER_HPP


#include "ProcessingClass.hpp"
#include <list>


namespace gpstk
{

      /** @addtogroup GPSsolutions */
      //@{


      /** This is a class to compute the mean and variance of the MW combination.
       *
       */
    class MWFilter : public ProcessingClass
    {
    public:

         /// Default constructor, setting default parameters.
        MWFilter()
        {};



         /** Return a satTypeValueMap object, adding the new data generated
          *  when calling this object.
          *
          * @param epoch     Time of observations.
          * @param gData     Data object holding the data.
          */
        virtual satTypeValueMap& Process( const CommonTime& epoch,
                                          satTypeValueMap& gData )
            throw(ProcessingException);


         /** Return a gnssSatTypeValue object, adding the new data generated
          *  when calling this object.
          *
          * @param gData    Data object holding the data.
          */
        virtual gnssSatTypeValue& Process(gnssSatTypeValue& gData)
            throw(ProcessingException)
        { (*this).Process(gData.header.epoch, gData.body); return gData; };


         /** Return a gnssRinex object, adding the new data generated when
          *  calling this object.
          *
          * @param gData    Data object holding the data.
          */
        virtual gnssRinex& Process(gnssRinex& gData)
            throw(ProcessingException);


         /** Return a gnssDataMap object, adding the new data generated when
          *  calling this object.
          *
          * @param gData    Data object holding the data.
          */
        virtual gnssDataMap& Process(gnssDataMap& gData)
            throw(ProcessingException);


         /// Return a string identifying this object.
        virtual std::string getClassName(void) const;


         /// Destructor
        virtual ~MWFilter() {};


    private:


        typedef std::map<SatID, double> ArcData;
        typedef std::map<SourceID, ArcData> ArcDataMap;

         /// Map holding the information regarding every satellite
        ArcData arcData;

         /// Map holding the information for gnssDataMap
        ArcDataMap arcDataMap;

         /// A structure used to store filter data for a SV.
        struct filterData
        {
            // Default constructor initializing the data in the structure
            filterData()
                : formerEpoch(CommonTime::BEGINNING_OF_TIME),
                  windowSize(0), meanMW(0.0), varMW(0.0)
            {};

            CommonTime formerEpoch; ///< The previous epoch time stamp.
            int windowSize;         ///< Size of current window, in samples.
            double meanMW;          ///< Accumulated mean value of combination.
            double varMW;           /// Variance
        };


        typedef std::map<SatID, filterData> MWData;
        typedef std::map<SourceID, MWData> MWDataMap;

         /// Map holding the information regarding every satellite
        MWData mwData;

         /// Map holding the information for gnssDataMap
        MWDataMap mwDataMap;


    }; // End of class 'MWFilter'

      //@}

}  // End of namespace gpstk

#endif   // GPSTK_MWFILTER_HPP
