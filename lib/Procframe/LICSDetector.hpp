#pragma ident "$Id$"

/**
 * @file LICSDetector.hpp
 * This is a class to detect cycle slips using LI observables.
 */

#ifndef LICSDETECTOR_HPP
#define LICSDETECTOR_HPP

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


#include "ProcessingClass.hpp"


namespace gpstk
{

      /** @addtogroup GPSsolutions */
      //@{


      /** This is a class to detect cycle slips using LI observables.
       *
       * This class is meant to be used with the GNSS data structures objects
       * found in "DataStructures" class.
       *
       * A typical way to use this class follows:
       *
       * @code
       *   RinexObsStream rin("ebre0300.02o");
       *
       *   gnssRinex gRin;
       *   ComputeLI getLI;
       *   LICSDetector markCSLI;
       *
       *   while(rin >> gRin)
       *   {
       *      gRin >> getLI >> markCSLI;
       *   }
       * @endcode
       *
       * The "LICSDetector" object will visit every satellite in the GNSS data
       * structure that is "gRin" and will decide if a cycle slip has happened
       * in the given observable.
       *
       * The algorithm will use LI observables, and the LLI1 and LLI2 indexes.
       * The result (a 1 if a cycle slip is found, 0 otherwise) will be stored
       * in the data structure both as the CSL1 and CSL2 indexes.
       *
       * This algorithm will use some values as maximum interval of time
       * between two successive epochs, minimum threshold to declare cycle
       * slip and LI combination limit drift.
       *
       * The default values are usually fine, but you may change them with the
       * appropriate methods. This is of special importance for the maximum
       * interval time, that should be adjusted for your sampling rate. It is
       * 61 seconds by default, which is appropriate for 30 seconds per sample
       * RINEX observation files.
       *
       * When used with the ">>" operator, this class returns the same incoming
       * data structure with the cycle slip indexes inserted along their
       * corresponding satellites. Be warned that if a given satellite does not
       * have the observations required, it will be summarily deleted from the
       * data structure.
       *
       * Be aware that some combinations of cycle slips in L1 and L2 may result
       * in a false negative when using a cycle slip detector based on LI.
       * Therefore, to be on the safe side you should complement this with
       * another kind of detector, such as one based on the Melbourne-Wubbena
       * combination.
       *
       * @sa MWCSDetector.hpp for more information.
       *
       * \warning Cycle slip detectors are objets that store their internal
       * state, so you MUST NOT use the SAME object to process DIFFERENT data
       * streams.
       *
       */
    class LICSDetector : public ProcessingClass
    {
    public:

         /// Default constructor, setting default parameters.
        LICSDetector()
            : checkExist(true),
              satSys(SatID::systemGPS), obsType(TypeID::LI12),
              lliType1(TypeID::LLI1), lliType2(TypeID::LLI2),
              resultType1(TypeID::CSL1), resultType2(TypeID::CSL2),
              deltaTMax(61.0), minThreshold(0.04), LIDrift(0.002),
              useLLI(true)
        {};


         /** Common constructor.
          *
          * @param mThr    Minimum threshold to declare cycle slip, in meters.
          * @param drift   LI combination limit drift, in meters/second.
          * @param dtMax   Maximum interval of time allowed between two
          *                successive epochs, in seconds.
          */
        LICSDetector( const double& mThr,
                      const double& drift,
                      const double& dtMax = 61.0,
                      const bool& use = true );


         /** Return a satTypeValueMap object, adding the new data generated
          *  when calling this object.
          *
          * @param epoch     Time of observations.
          * @param gData     Data object holding the data.
          * @param epochflag Epoch flag.
          */
        virtual satTypeValueMap& Process( const CommonTime& epoch,
                                          satTypeValueMap& gData,
                                          const short& epochflag = 0 )
            throw(ProcessingException);


        /** Method to set the satSys.
         *
         * @param sys      SatSys
         */
        virtual LICSDetector& setSatSystem(const SatID::SatelliteSystem& sys)
        { satSys = sys; return (*this); };


        /** Method to get the satSys.
         */
        virtual SatID::SatelliteSystem getSatSystem() const
        { return satSys; };


         /** Method to get the obsType.
          */
        virtual TypeID getObsType() const
        { return obsType; };


         /** Method to set the obsType.
          *
          * @param type      ObsType
          */
        virtual LICSDetector& setObsType(const TypeID& type)
        { obsType = type; return (*this); };


         /** Method to get the lliType1.
          */
        virtual TypeID getLLIType1() const
        { return lliType1; };


         /** Method to set the lliType1.
          *
          * @param type      LLIType1
          */
        virtual LICSDetector& setLLIType1(const TypeID& type)
        { lliType1 = type; return (*this); };


         /** Method to get the lliType2.
          */
        virtual TypeID getLLIType2() const
        { return lliType2; };


         /** Method to set the lliType2.
          *
          * @param type      LLIType2
          */
        virtual LICSDetector& setLLIType2(const TypeID& type)
        { lliType2 = type; return (*this); };


         /** Method to get the resultType1.
          */
        virtual TypeID getResultType1() const
        { return resultType1; };


         /** Method to set the resultType1.
          *
          * @param type      ResultType1
          */
        virtual LICSDetector& setResultType1(const TypeID& type)
        { resultType1 = type; return (*this); };


         /** Method to get the resultType2.
          */
        virtual TypeID getResultType2() const
        { return resultType2; };


         /** Method to set the resultType2.
          *
          * @param type      ResultType2
          */
        virtual LICSDetector& setResultType2(const TypeID& type)
        { resultType2 = type; return (*this); };


         /** Method to get the maximum interval of time allowed between two
          *  successive epochs, in seconds.
          */
        virtual double getDeltaTMax() const
        { return deltaTMax; };


         /** Method to set the maximum interval of time allowed between two
          *  successive epochs.
          *
          * @param maxDelta      Maximum interval of time, in seconds
          */
        virtual LICSDetector& setDeltaTMax(const double& maxDelta)
        {
            if (maxDelta > 0.0)
            {
                deltaTMax = maxDelta;
            }
            else
            {
                deltaTMax = 61.0;
            }

            return (*this);
        };


         /** Method to get the minimum threshold for cycle slip detection, in
          *  meters.
          */
        virtual double getMinThreshold() const
        { return minThreshold; };


         /** Method to set the minimum threshold for cycle slip detection, in
          *  meters.
          *
          * @param mThr    Minimum threshold for cycle slip detection, in
          *                meters.
          */
        virtual LICSDetector& setMinThreshold(const double& mThr)
        {
            if (mThr < 0.0)
            {
                minThreshold = 0.04;
            }
            else
            {
                minThreshold = mThr;
            }

            return (*this);
        };


         /// Method to get the LI combination limit drift, in meters/second
        virtual double getLIDrift() const
        { return LIDrift; };


         /** Method to set the LI combination limit drift, in meters/second
          *
          * @param drift     LI combination limit drift, in meters/second.
          */
        virtual LICSDetector& setLIDrift(const double& drift)
        {
            if (drift > 0.0)
            {
                LIDrift = drift;
            }
            else
            {
                LIDrift = 0.002;
            }

            return (*this);
        };


         /// Method to know if the LLI check is enabled or disabled.
        virtual bool getUseLLI() const
        { return useLLI; };


         /** Method to set whether the LLI indexes will be used as an aid
          *  or not.
          *
          * @param use   Boolean value enabling/disabling LLI check
          */
        virtual LICSDetector& setUseLLI(const bool& use)
        { useLLI = use; return (*this); };


         /** Return a gnssSatTypeValue object, adding the new data generated
          *  when calling this object.
          *
          * @param gData    Data object holding the data.
          */
        virtual gnssSatTypeValue& Process(gnssSatTypeValue& gData)
            throw(ProcessingException)
        { Process(gData.header.epoch, gData.body); return gData; };


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


        virtual LICSDetector& setCheckExist(const bool& check)
        { checkExist = check; return (*this); };


        virtual bool getCheckExist() const
        { return checkExist; };


         /// Return a string identifying this object.
        virtual std::string getClassName(void) const;


         /// Destructor
        virtual ~LICSDetector() {};


    private:

        bool checkExist;

        /// System of satellite.
        SatID::SatelliteSystem satSys;

        /// Type of observable.
        TypeID obsType;

        /// Type of LLI1 record.
        TypeID lliType1;

        /// Type of LLI2 record.
        TypeID lliType2;

        /// Type of result #1.
        TypeID resultType1;

        /// Type of result #2.
        TypeID resultType2;


        /// Maximum interval of time allowed between two successive epochs,
        /// in seconds.
        double deltaTMax;

        /// Minimum threshold to declare cycle slip, in meters.
        double minThreshold;

        /// LI combination limit drift, in meters/second.
        double LIDrift;

        /// This field tells whether to use or ignore the LLI indexes as
        /// an aid.
        bool useLLI;


        /// A structure used to store filter data for a SV.
        struct filterData
        {
            // Default constructor initializing the data in the structure
            filterData() : formerEpoch(CommonTime::BEGINNING_OF_TIME),
                           windowSize(0), formerLI(0.0), formerBias(0.0),
                           formerDeltaT(1.0)
            {};

            CommonTime formerEpoch; ///< The previous epoch time stamp.
            int windowSize;         ///< Size of current window, in samples.
            double formerLI;        ///< Value of the previous LI observable.
            double formerBias;      ///< Previous bias (LI_1 - LI_0).
            double formerDeltaT;    ///< Previous time difference, in seconds.
        };


        typedef std::map<SatID, filterData> LIData;
        typedef std::map<SourceID, LIData> LIDataMap;

        /// Map holding the information regarding every satellite
        LIData liData;

        /// Map holding the information for gnssDataMap
        LIDataMap liDataMap;


         /** Method that implements the LI cycle slip detection algorithm.
          *
          * @param epoch     Time of observations.
          * @param sat       SatID.
          * @param tvMap     Data structure of TypeID and values.
          * @param epochflag Epoch flag.
          * @param li        Current LI observation value.
          * @param lli1      LLI1 index.
          * @param lli2      LLI2 index.
          */
        virtual double getDetection( const CommonTime& epoch,
                                     const SatID& sat,
                                     typeValueMap& tvMap,
                                     const short& epochflag,
                                     const double& li,
                                     const double& lli1,
                                     const double& lli2 );

    }; // End of class 'LICSDetector'

    //@}

}  // End of namespace gpstk

#endif   // LICSDETECTOR_HPP
