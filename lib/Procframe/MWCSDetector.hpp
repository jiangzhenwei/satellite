#pragma ident "$Id$"

/**
 * @file MWCSDetector.hpp
 * This is a class to detect cycle slips using the Melbourne-Wubbena
 * combination.
 */

#ifndef MWCSDETECTOR_HPP
#define MWCSDETECTOR_HPP

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
//  Hui Wei, 2016
//
//============================================================================



#include "ProcessingClass.hpp"
#include <list>


namespace gpstk
{

      /** @addtogroup GPSsolutions */
      //@{


      /** This is a class to detect cycle slips using MW observables.
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
       *   ComputeMelbourneWubbena getMW;
       *   MWCSDetector markCSMW;
       *
       *   while(rin >> gRin)
       *   {
       *      gRin >> getMW >> markCSMW;
       *   }
       * @endcode
       *
       * The "MWCSDetector" object will visit every satellite in the GNSS data
       * structure that is "gRin" and will decide if a cycle slip has happened
       * in the given observable.
       *
       * The algorithm will use MW observables, and the LLI1 and LLI2 indexes.
       * The result (a 1 if a cycle slip is found, 0 otherwise) will be stored
       * in the data structure both as the CSL1 and CSL2 indexes.
       *
       * In taking the decision, this algorithm will use criteria as the
       * maximum interval of time between two successive epochs and the
       * maximum number of Melbourne-Wubbena wavelenghts allowed above or
       * below the MW combination average for that arc.
       *
       * The default values are usually fine, but you may change them with the
       * appropriate methods. This is of special importance for the maximum
       * interval time, that should be adjusted for your sampling rate. It is
       * 61 seconds by default, which is appropriate for 30 seconds per sample
       * RINEX observation files.
       *
       * When used with the ">>" operator, this class returns the same
       * incoming data structure with the cycle slip indexes inserted along
       * their corresponding satellites. Be warned that if a given satellite
       * does not have the observations required, it will be summarily deleted
       * from the data structure.
       *
       * You should be aware that the Melbourne-Wubbena combination is based
       * on a mix of code and phase observations, so it is very noisy.
       * Therefore, it has a tendency to yield a high number of false
       * positives if you are not careful with its parameters. Because of
       * this, the default parameters are very conservative, i.e., the
       * detector is NOT very sensitive by default.
       *
       * Best results are achieved when using this detector as a "backup"
       * detector for detectors based in LI combination, like this:
       *
       * @code
       *   RinexObsStream rin("ebre0300.02o");
       *
       *   gnssRinex gRin;
       *   ComputeLI getLI;
       *   LICSDetector markCSLI;
       *   ComputeMelbourneWubbena getMW;
       *   MWCSDetector markCSMW;
       *
       *   while(rin >> gRin)
       *   {
       *      gRin >> getLI >> getMW >> markCSLI >> markCSMW;
       *   }
       * @endcode
       *
       * @sa LICSDetector.hpp and LICSDetector2.hpp for more information.
       *
       * \warning Cycle slip detectors are objets that store their internal
       * state, so you MUST NOT use the SAME object to process DIFFERENT data
       * streams.
       *
       */
    class MWCSDetector : public ProcessingClass
    {
    public:

        /// Default constructor, setting default parameters.
        MWCSDetector()
            : checkExist(true),
              satSys(SatID::systemGPS), obsType(TypeID::MW12),
              lliType1(TypeID::LLI1), lliType2(TypeID::LLI2),
              resultType1(TypeID::CSL1), resultType2(TypeID::CSL2),
              deltaTMax(61.0), maxNumLambdas(10.0*WL_WAVELENGTH_GPS_L1L2),
              useLLI(true)
        {};


        /** Common constructor.
         *
         * @param mLambdas      Maximum deviation allowed before declaring
         *                      cycle slip (in number of Melbourne-Wubbena
         *                      wavelenghts).
         * @param dtMax         Maximum interval of time allowed between two
         *                      successive epochs, in seconds.
         */
        MWCSDetector( const double& mLambdas,
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
        virtual MWCSDetector& setSatSystem(const SatID::SatelliteSystem& sys)
        { satSys = sys; return (*this); };


        /** Method to get the satSys.
         */
        virtual SatID::SatelliteSystem getSatSystem() const
        { return satSys; };


        /** Method to set the obsType.
         *
         * @param type      ObsType
         */
        virtual MWCSDetector& setObsType(const TypeID& type)
        { obsType = type; return (*this); };


        /** Method to get the obsType.
         */
        virtual TypeID getObsType() const
        { return obsType; };


        /** Method to set the lliType1.
         *
         * @param type      lliType1
         */
        virtual MWCSDetector& setLLIType1(const TypeID& type)
        { lliType1 = type; return (*this); };


        /** Method to get the lliType1.
         */
        virtual TypeID getLLIType1() const
        { return lliType1; };


        /** Method to set the lliType2.
         *
         * @param type      lliType2
         */
        virtual MWCSDetector& setLLIType2(const TypeID& type)
        { lliType2 = type; return (*this); };


        /** Method to get the lliType2.
         */
        virtual TypeID getLLIType2() const
        { return lliType2; };


        /** Method to set the resultType1.
         *
         * @param type      resultType1
         */
        virtual MWCSDetector& setResultType1(const TypeID& type)
        { resultType1 = type; return (*this); };


        /** Method to get the resultType1.
         */
        virtual TypeID getResultType1() const
        { return resultType1; };


        /** Method to set the resultType2.
         *
         * @param type      resultType2
         */
        virtual MWCSDetector& setResultType2(const TypeID& type)
        { resultType2 = type; return (*this); };


        /** Method to get the resultType2.
         */
        virtual TypeID getResultType2() const
        { return resultType2; };


        /** Method to set the maximum interval of time allowed between two
         *  successive epochs.
         *
         * @param maxDelta      Maximum interval of time, in seconds
         */
        virtual MWCSDetector& setDeltaTMax(const double& maxDelta)
        {
            if(maxDelta > 0.0)
            {
                deltaTMax = maxDelta;
            }
            else
            {
                deltaTMax = 61.0;
            }

            return (*this);
        };


        /** Method to get the maximum interval of time allowed between two
         *  successive epochs, in seconds.
         */
        virtual double getDeltaTMax() const
        { return deltaTMax; };


        /** Method to set the maximum deviation allowed before declaring
         *  cycle slip (in number of Melbourne-Wubbena wavelenghts).
         *
         * @param mLambdas     Maximum deviation allowed before declaring
         *                     cycle slip (in number of Melbourne-Wubbena
         *                     wavelenghts).
         */
        virtual MWCSDetector& setMaxNumLambdas(const double& mLambdas)
        {
            if(mLambdas > 0.0)
            {
                maxNumLambdas = mLambdas;
            }
            else
            {
                maxNumLambdas = 10.0*WL_WAVELENGTH_GPS_L1L2;
            }

            return (*this);
        };


        /** Method to get the maximum deviation allowed before declaring
         *  cycle slip (in number of Melbourne-Wubbena wavelenghts).
         */
        virtual double getMaxNumLambdas() const
        { return maxNumLambdas; };


        /** Method to set whether the LLI indexes will be used as
         *  an aid or not.
         *
         * @param use   Boolean value enabling/disabling LLI check.
         */
        virtual MWCSDetector& setUseLLI(const bool& use)
        { useLLI = use; return (*this); };


        /// Method to know if the LLI check is enabled or disabled.
        virtual bool getUseLLI() const
        { return useLLI; };


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


        virtual MWCSDetector& setCheckExist(const bool& check)
        { checkExist = check; return (*this); };


        virtual bool getCheckExist() const
        { return checkExist; };


        /// Return a string identifying this object.
        virtual std::string getClassName(void) const;


        /// Destructor
        virtual ~MWCSDetector() {};


    private:

        bool checkExist;

        /// System of satellite.
        SatID::SatelliteSystem satSys;

        /// Type of observation.
        TypeID obsType;

        /// Type of LMW1 record.
        TypeID lliType1;

        /// Type of LMW2 record.
        TypeID lliType2;

        /// Type of result #1.
        TypeID resultType1;

        /// Type of result #2.
        TypeID resultType2;


        /** Maximum interval of time allowed between two successive
         *  epochs, in seconds.
         */
        double deltaTMax;

        /** Maximum deviation allowed before declaring cycle slip,
         *  in number of Melbourne-Wubbena wavelenghts.
         */
        double maxNumLambdas;

        /// Whether use or ignore the LLI indexes as an aid.
        bool useLLI;


        /// A structure used to store filter data for a SV.
        struct filterData
        {
            // Default constructor initializing the data in the structure
            filterData() : formerEpoch(CommonTime::BEGINNING_OF_TIME),
                           windowSize(0), meanMW(0.0)
            {};

            CommonTime formerEpoch; ///< The previous epoch time stamp.
            int windowSize;         ///< Size of current window, in samples.
            double meanMW;          ///< Accumulated mean value of combination.
            double varMW;           ///< Accumulated std value of combination.
        };

        typedef std::map<SatID, filterData> MWData;
        typedef std::map<SourceID, MWData> MWDataMap;

        /// Map holding the information regarding every satellite
        MWData mwData;

        /// Map holding the information for gnssDataMap
        MWDataMap mwDataMap;


        /** Method that implements the Melbourne-Wubbena cycle slip
         *  detection algorithm.
         *
         * @param epoch     Time of observations.
         * @param sat       SatID.
         * @param tvMap     Data structure of TypeID and values.
         * @param epochflag Epoch flag.
         * @param mw        Current MW observation value.
         * @param lli1      LLI1 index.
         * @param lli2      LLI2 index.
         */
        virtual double getDetection( const CommonTime& epoch,
                                     const SatID& sat,
                                     typeValueMap& tvMap,
                                     const short& epochflag,
                                     const double& mw,
                                     const double& lli1,
                                     const double& lli2 );


    }; // End of class 'MWCSDetector'

    //@}

}  // End of namespace gpstk

#endif   // MWCSDETECTOR_HPP
