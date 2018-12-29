#pragma ident "$Id$"

/**
 * @file ComputePartials.hpp
 */

#ifndef BASIC_MODEL2_HPP
#define BASIC_MODEL2_HPP

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



#include "ProcessingClass.hpp"
#include "ReferenceSystem.hpp"
#include "MSCStore.hpp"


namespace gpstk
{
      /** @addtogroup GPSsolutions */
      //@{

    class ComputePartials : public ProcessingClass
    {
    public:

        /// Default constructor.
        ComputePartials()
            : pRefSys(NULL), pMSCStore(NULL)
        {};


        /// Explicit constructor.
        ComputePartials( const satVectorMap& state,
                         ReferenceSystem& refSys,
                         MSCStore& mscStore )
        {
            satStateMap = state;
            pRefSys = &refSys;
            pMSCStore = &mscStore;
        };


         /** Return a satTypeValueMap object, adding the new data generated
          *  when calling a modeling object.
          *
          * @param time      Epoch.
          * @param gData     Data object holding the data.
          */
        virtual satTypeValueMap& Process(const CommonTime& time,
                                         satTypeValueMap& gData)
            throw(ProcessingException);


         /** Return a gnssSatTypeValue object, adding the new data generated
          *  when calling a modeling object.
          *
          * @param gData    Data object holding the data.
          */
        virtual gnssSatTypeValue& Process(gnssSatTypeValue& gData)
            throw(ProcessingException)
        { Process(gData.header.epoch, gData.body); return gData; };


         /** Return a gnssRinex object, adding the new data generated when
          *  calling a modeling object.
          *
          * @param gData    Data object holding the data.
          */
        virtual gnssRinex& Process(gnssRinex& gData)
            throw(ProcessingException)
        { Process(gData.header.epoch, gData.body); return gData; };


         /** Return a gnssDataMap object, adding the new data generated when
          *  calling a modeling object.
          *
          * @param gData    Data object holding the data.
          */
        virtual gnssDataMap& Process(gnssDataMap& gData)
            throw(ProcessingException);


         /// Method to get sat state to be used.
        virtual satVectorMap getSatStateMap() const
        { return satStateMap; };


         /** Method to set the sat orbit to be used.
          *
          * @param state    satVectorMap object to be used
          */
        virtual ComputePartials& setSatStateMap(const satVectorMap& state)
        { satStateMap = state; return (*this); };


         /// Method to get source clock to be used.
        virtual sourceValueMap getSourceClockMap() const
        { return sourceClockMap; };


         /** Method to set the station clock to be used.
          *
          * @param clockStore   sourceValueMap object to be used
          */
        virtual ComputePartials& setSourceClockMap(const sourceValueMap& clockStore)
        { sourceClockMap = clockStore; return (*this); };


         /// Method to get a pointer to the ReferenceSystem to be used
         /// with GNSS data structures.
        virtual ReferenceSystem* getReferenceSystem() const
        { return pRefSys; };


         /** Method to set the ReferenceSystem to be used with GNSS
          *  data structures.
          *
          * @param refSys     ReferenceSystem object to be used
          */
        virtual ComputePartials& setReferenceSystem(ReferenceSystem& refSys)
        { pRefSys = &refSys; return (*this); };


         /// Return a pointer to the MSCStore object currently in use.
        virtual MSCStore *getMSCStore(void) const
        { return pMSCStore; };


         /** Sets MSCStore object to be used.
          *
          * @param msc  MSCStore object.
          */
        virtual ComputePartials& setMSCStore(MSCStore& msc)
        { pMSCStore = &msc; return (*this); };


        /// Get satellite pos in ECI.
        virtual satVectorMap getSatPosECIMap() const
        { return satPosECIMap; };


        /// Get satellite pos in ECF.
        virtual satVectorMap getSatPosECFMap() const
        { return satPosECFMap; };


        /// Get source pos in ECI.
        virtual sourceVectorMap getSourcePosECIMap() const
        { return sourcePosECIMap; };


        /// Get source pos in ECF.
        virtual sourceVectorMap getSourcePosECFMap() const
        { return sourcePosECFMap; };


         /// Return a string identifying this object.
        virtual std::string getClassName(void) const;


         /// Destructor.
        virtual ~ComputePartials() {};


    protected:

        /// Satellite state
        satVectorMap satStateMap;

        satVectorMap satPosECIMap;
        satVectorMap satPosECFMap;

        Vector<double> posSatECI;
        Vector<double> posSatECF;

        Vector<double> velSatECI;
        Vector<double> velSatECF;

        /// Source Clock
        sourceValueMap sourceClockMap;

        double cdtSource;

        sourceVectorMap sourcePosECIMap;
        sourceVectorMap sourcePosECFMap;

        Vector<double> posSourceECI;
        Vector<double> posSourceECF;

        /// ReferenceSystem
        ReferenceSystem* pRefSys;

        /// Source Position
        MSCStore* pMSCStore;

        Position nominalPos;

    }; // End of class 'ComputePartials'

    //@}

}  // End of namespace gpstk

#endif   // COMPUTE_PARTIALS_HPP

