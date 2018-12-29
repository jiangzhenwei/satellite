#pragma ident "$Id$"

/**
 * @file ComputeWindUp.hpp
 * This class computes the wind-up effect on the phase observables, in radians.
 */

#ifndef COMPUTE_WINDUP_HPP
#define COMPUTE_WINDUP_HPP

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


#include <string>
#include "ProcessingClass.hpp"
#include "ReferenceSystem.hpp"
#include "SolarSystem.hpp"
#include "MSCStore.hpp"
#include "ComputeSatAttitude.hpp"
#include "constants.hpp"


namespace gpstk
{

      /** @addtogroup DataStructures */
      //@{


      /** This class computes the wind-up effect on the phase observables,
       *  in radians.
       *
       * This class is meant to be used with the GNSS data structures objects
       * found in "DataStructures" class.
       *
       * A typical way to use this class follows:
       *
       * @code
       * @endcode
       *
       * The "ComputeWindUp" object will visit every satellite in the GNSS
       * data structure that is "gRin" and will compute the corresponding
       * receiver-satellite wind-up effect, in radians.
       *
       * When used with the ">>" operator, this class returns the same
       * incoming data structure with the wind-up inserted in it. Be warned
       * that if a given satellite does not have the observations required,
       * it will be summarily deleted from the data structure.
       *
       * \warning ComputeWindUp objects store their internal state, so
       * you MUST NOT use the SAME object to process DIFFERENT data streams.
       *
       * \warning It is recommended that ComputeWindUp objects are used after
       * calling a SatArcMarker object, because they work better when cycle
       * slips are properly tracked.
       *
       */
    class ComputeWindUp : public ProcessingClass
    {
    public:

         /// Default constructor
        ComputeWindUp()
            : pRefSys(NULL), pSolSys(NULL),
              pEphStore(NULL), pMSCStore(NULL)
        {
            nominalPos = Position(0.0,0.0,0.0,Position::Cartesian,NULL);
        };


         /** Return a satTypeValueMap object, adding the new data generated
          *  when calling this object.
          *
          * @param time      Epoch corresponding to the data.
          * @param gData     Data object holding the data.
          */
        virtual satTypeValueMap& Process( const CommonTime& time,
                                          satTypeValueMap& gData )
            throw(ProcessingException);


         /** Return a gnssSatTypeValue object, adding the new data
          *  generated when calling this object.
          *
          * @param gData    Data object holding the data.
          */
        virtual gnssSatTypeValue& Process(gnssSatTypeValue& gData)
            throw(ProcessingException)
        { Process(gData.header.epoch, gData.body); return gData; };


         /** Return a gnssRinex object, adding the new data generated
          *  when calling this object.
          *
          * @param gData    Data object holding the data.
          */
        virtual gnssRinex& Process(gnssRinex& gData)
            throw(ProcessingException)
        { Process(gData.header.epoch, gData.body); return gData; };


         /** Return a gnssDataMap object, adding the new data generated
          *  when calling this object.
          *
          * @param gData    Data object holding the data.
          */
        virtual gnssDataMap& Process(gnssDataMap& gData)
            throw(ProcessingException);


        /** Set ReferenceSystem object to be used.
        *
        * @param refSys  ReferenceSystem object.
        */
        virtual ComputeWindUp& setReferenceSystem(ReferenceSystem& refSys)
        { pRefSys = &refSys; return (*this); };


        /** Set SolarSystem object to be used.
        *
        * @param solSys  SolarSystem object.
        */
        virtual ComputeWindUp& setSolarSystem(SolarSystem& solSys)
        { pSolSys = &solSys; return (*this); };


        /** Set ephemeris object to be used.
         *
         * @param ephStore      Ephemeris object.
         */
        virtual ComputeWindUp& setEphStore(XvtStore<SatID>& ephStore)
        { pEphStore = &ephStore; return (*this); };


         /// Return nominal position of station.
        virtual Position getNominalPosition(void) const
        { return nominalPos; };


         /** Set nominal position of station.
          * @param staPos    Nominal position of station.
          */
        virtual ComputeWindUp& setNominalPosition(const Position& staPos)
        { nominalPos = staPos; return (*this); };


         /// Return a pointer to the MSCStore object currently in use.
        virtual MSCStore *getMSCStore(void) const
        { return pMSCStore; };


         /** Sets MSCStore object to be used.
          *
          * @param msc  MSCStore object.
          */
        virtual ComputeWindUp& setMSCStore(MSCStore& mscStore)
        { pMSCStore = &mscStore; return (*this); };


         /** Sets satYawDataMap to be used.
          *
          * @param yawData  satYawDataMap.
          */
        virtual ComputeWindUp& setSatYawData(const satYawDataMap& yawData)
        { satYawData = yawData; return (*this); };


         /// Return a string identifying this object.
        virtual std::string getClassName(void) const;


         /// Destructor
        virtual ~ComputeWindUp() {};


    private:

        /// Pointer to ReferenceSystem
        ReferenceSystem* pRefSys;

        /// Pointer to SolarSystem
        SolarSystem* pSolSys;

        /// Pointer to ephemeris object
        XvtStore<SatID> *pEphStore;

        /// Pointer to MSCStore object
        MSCStore* pMSCStore;

        /// Station position
        Position nominalPos;

        /// Sat yaw data
        satYawDataMap satYawData;


        /// A structure used to store phase data.
        struct phaseData
        {
            // Default constructor initializing the data in the structure
            phaseData() : satPreviousPhase(0.0),
                          staPreviousPhase(0.0),
                          arcNum(0.0)
            {};

            double satPreviousPhase;      ///< Previous phase for Sat.
            double staPreviousPhase;      ///< Previous phase for Station.
            double arcNum;                ///< Satellite arc number.
        };


        typedef std::map<SatID, phaseData> SatPhaseData;
        typedef std::map<SourceID, SatPhaseData> SatPhaseDataMap;

        /// Map to store satellite phase data
        SatPhaseData satPhaseData;

        /// Map to store station phase date
        SatPhaseDataMap satPhaseDataMap;

    }; // End of class 'ComputeWindUp'

      //@}

}  // End of namespace gpstk

#endif // COMPUTE_WINDUP_HPP
