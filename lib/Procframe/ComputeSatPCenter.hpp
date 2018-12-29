#pragma ident "$Id$"

/**
 * @file ComputeSatPCenter.hpp
 * This class computes the satellite antenna phase correction, in meters.
 */

#ifndef COMPUTE_SATPCENTER_HPP
#define COMPUTE_SATPCENTER_HPP

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


#include <cmath>
#include <string>
#include <sstream>
#include "ProcessingClass.hpp"
#include "ReferenceSystem.hpp"
#include "SolarSystem.hpp"
#include "SP3EphemerisStore.hpp"
#include "MSCStore.hpp"
#include "AntexReader.hpp"
#include "ComputeSatAttitude.hpp"
#include "StringUtils.hpp"
#include "constants.hpp"


namespace gpstk
{

      /** @addtogroup DataStructures */
      //@{


      /** This class computes the satellite antenna phase correction, in meters.
       *
       * This class is meant to be used with the GNSS data structures objects
       * found in "DataStructures" class.
       *
       * A typical way to use this class follows:
       *
       * @code
       * @endcode
       *
       * The "ComputeSatPCenter" object will visit every satellite in the GNSS
       * data structure that is "gRin" and will compute the corresponding
       * satellite antenna phase correction, in meters.
       *
       * When used with the ">>" operator, this class returns the same
       * incoming data structure with the "satPCenter" TypeID inserted in it.
       * Be warned that if a given satellite does not have the required data,
       * it will be summarily deleted from the data structure.
       *
       * \warning The ComputeSatPCenter objects generate corrections that are
       * interpreted as an "advance" in the signal, instead of a delay.
       * Therefore, those corrections always hava a negative sign.
       *
       */
    class ComputeSatPCenter : public ProcessingClass
    {
    public:

         /// Default constructor
        ComputeSatPCenter()
            : pRefSys(NULL), pSolSys(NULL),
              pEphStore(NULL), pMSCStore(NULL),
              pAntexReader(NULL)
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


         /// Return a pointer to the AntexReader object currently in use.
        virtual AntexReader *getAntexReader(void) const
        { return pAntexReader; };


         /** Sets AntexReader object to be used.
          *
          * @param antexObj  AntexReader object containing satellite
          *                  antenna data.
          */
        virtual ComputeSatPCenter& setAntexReader(AntexReader& antexReader)
        { pAntexReader = &antexReader; return (*this); };


        /** Set ReferenceSystem object to be used.
        *
        * @param ref  ReferenceSystem object.
        */
        virtual ComputeSatPCenter& setReferenceSystem(ReferenceSystem& ref)
        { pRefSys = &ref; return (*this); };


        /** Set SolarSystem object to be used.
        *
        * @param sol  SolarSystem object.
        */
        virtual ComputeSatPCenter& setSolarSystem(SolarSystem& sol)
        { pSolSys = &sol; return (*this); };


        /// Return a pointer to the ephemeris object currently in use.
        virtual XvtStore<SatID>* getEphStore(void) const
        { return pEphStore; };


        /** Set ephemeris object to be used.
         *
         * @param ephStore      Ephemeris object.
         */
        virtual ComputeSatPCenter& setEphStore(XvtStore<SatID>& ephStore)
        { pEphStore = &ephStore; return (*this); };


         /// Return nominal position of station.
        virtual Position getNominalPosition(void) const
        { return nominalPos; };


         /** Set nominal position of station.
          * @param staPos    Nominal position of station.
          */
        virtual ComputeSatPCenter& setNominalPosition(const Position& staPos)
        { nominalPos = staPos; return (*this); };


        /// Return a pointer to the MSCStore object currently in use.
        virtual MSCStore* getMSCStore(void) const
        { return pMSCStore; };


         /** Set MSCStore object to be used.
          *
          * @param mscStore     MSCStore object.
          */
        virtual ComputeSatPCenter& setMSCStore(MSCStore& mscStore)
        { pMSCStore = &mscStore; return (*this); };


         /** Set satYawDataMap to be used.
          *
          * @param yawData     satYawDataMap.
          */
        virtual ComputeSatPCenter& setSatYawData(const satYawDataMap& yawData)
        { satYawData = yawData; return (*this); };


         /// Return a string identifying this object.
        virtual std::string getClassName(void) const;


         /// Destructor
        virtual ~ComputeSatPCenter() {};


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

        /// Pointer to AntexReader object
        AntexReader* pAntexReader;

        /// Sat yaw
        satYawDataMap satYawData;

    }; // End of class 'ComputeSatPCenter'

      //@}

}  // End of namespace gpstk

#endif // COMPUTE_SATPCENTER_HPP
