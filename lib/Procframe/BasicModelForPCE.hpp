#pragma ident "$Id$"

/**
 * @file BasicModelForPCE.hpp
 * This is a class to compute the basic parts of a GNSS model, i.e.:
 * Geometric distance, relativity correction, satellite position and
 * velocity at transmission time, satellite elevation and azimuth, etc.
 */

#ifndef BASICMODEL_FOR_PCE_HPP
#define BASICMODEL_FOR_PCE_HPP

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
#include "EphemerisRange.hpp"
#include "XvtStore.hpp"
#include "MSCStore.hpp"


namespace gpstk
{
      /** @addtogroup GPSsolutions */
      //@{

      /** This is a class to compute the basic parts of a GNSS model, like
       *  geometric distance, relativity correction, satellite position and
       *  velocity at transmission time, satellite elevation and azimuth, etc.
       *
       * This class is intended to be used with GNSS Data Structures (GDS).
       *
       * A typical way to use this class follows:
       *
       * @code
       *      // Input observation file stream
       *   RinexObsStream rin("ebre0300.02o");
       *
       *      // Load the precise ephemeris file
       *   SP3EphemerisStore sp3Eph;
       *   sp3Eph.loadFile("igs11513.sp3");
       *
       *      // Reference position of receiver station
       *   Position nominalPos(4833520.2269, 41537.00768, 4147461.489);
       *
       *      // Some more code and definitions here...
       *
       *   gnssRinex gRin;  // GNSS data structure for fixed station data
       *
       *      // Set defaults of models. A typical C1-based modeling is used
       *   BasicModelForPCE model( nominalPos, sp3Eph );
       *
       *   while(rin >> gRin)
       *   {
       *
       *         // Apply the model on the GDS
       *      gRin >> model;
       *   }
       *
       * @endcode
       *
       * The "BasicModelForPCE" object will visit every satellite in
       * the GNSS data structure that is "gRin" and will try to compute
       * its model: Geometric distance, relativity delay, satellite position
       * at transmission time, satellite elevation and azimuth, etc.
       *
       * When used with the ">>" operator, this class returns the same
       * incoming data structure with the extra data inserted along their
       * corresponding satellites. Be warned that if a given satellite does
       * not have ephemeris information, it will be summarily deleted
       * from the data structure.
       *
       */
    class BasicModelForPCE : public ProcessingClass
    {
    public:

         /// Default constructor. Observable C1 will be used for computations
         /// and satellites with elevation less than 10 degrees will be
         /// deleted.
        BasicModelForPCE()
            : minElev(10.0),
              pSP3Store(NULL), pBCEStore(NULL), pMSCStore(NULL),
              defaultObsOfGPS(TypeID::C1), useTGDOfGPS(false),
              defaultObsOfGAL(TypeID::C1), useTGDOfGAL(false),
              defaultObsOfBDS(TypeID::C2), useTGDOfBDS(false)
        {
            nominalPos = Position(0.0,0.0,0.0,Position::Cartesian,NULL);
        };


         /** Explicit constructor taking as input reference
          *  station coordinates.
          *
          * Those coordinates may be Cartesian (X, Y, Z in meters) or Geodetic
          * (Latitude, Longitude, Altitude), but defaults to Cartesian.
          *
          * Also, a pointer to EllipsoidModel may be specified, but default is
          * NULL (in which case WGS84 values will be used).
          *
          * @param aRx   first coordinate [ X(m), or latitude (degrees N) ]
          * @param bRx   second coordinate [ Y(m), or longitude (degrees E) ]
          * @param cRx   third coordinate [ Z, height above ellipsoid or
          *              radius, in meters ]
          * @param s     coordinate system (default is Cartesian, may be set
          *              to Geodetic).
          * @param ell   pointer to EllipsoidModel.
          * @param frame Reference frame associated with this position.
          */
        BasicModelForPCE( const double& aRx,
                          const double& bRx,
                          const double& cRx,
                          Position::CoordinateSystem s = Position::Cartesian,
                          EllipsoidModel *ell = NULL,
                          ReferenceFrame frame = ReferenceFrame::Unknown );


        /// Explicit constructor, taking as input a Position object
        /// containing reference station coordinates.
        BasicModelForPCE(const Position& staCoordinates);


        /** Explicit constructor, taking as input reference station
         *  coordinates, ephemeris to be used and whether TGD will
         *  be computed or not.
         *
         * @param RxCoordinates Reference station coordinates.
         * @param dEphemeris    EphemerisStore object to be used by default.
         * @param dObservable   Observable type to be used by default.
         * @param applyTGD      Whether or not C1 observable will be
         *                      corrected from TGD effect.
         *
         */
        BasicModelForPCE( const Position& StaCoordinates,
                          XvtStore<SatID>& sp3Store,
                          XvtStore<SatID>& bceStore,
                          MSCStore& mscStore,
                          const TypeID& dObsOfGPS = TypeID::C1,
                          const TypeID& dObsOfGAL = TypeID::C1,
                          const TypeID& dObsOfBDS = TypeID::C2,
                          const bool& applyTGDOfGPS = false,
                          const bool& applyTGDOfGAL = false,
                          const bool& applyTGDOfBDS = false );


        /** Return a satTypeValueMap object, adding the new data generated
         *  when calling a modeling object.
         *
         * @param time      Epoch.
         * @param gData     Data object holding the data.
         */
        virtual satTypeValueMap& Process( const CommonTime& time,
                                          satTypeValueMap& gData )
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


        /// Method to get satellite elevation cut-off angle. By default, it
        /// is set to 10 degrees.
        virtual double getMinElev() const
        { return minElev; };


        /// Method to set satellite elevation cut-off angle. By default, it
        /// is set to 10 degrees.
        virtual BasicModelForPCE& setMinElev(double newElevation)
        { minElev = newElevation; return (*this); };


        /// Method to get the default observable for computations.
        virtual TypeID getDefaultObs( const SatID::SatelliteSystem& sys
                                                    = SatID::systemGPS ) const
        {
            if(sys == SatID::systemGPS)
                return defaultObsOfGPS;
            else if(sys == SatID::systemGalileo)
                return defaultObsOfGAL;
            else if(sys == SatID::systemBDS)
                return defaultObsOfBDS;
        };


        /** Method to set the default observable for computations.
         *
         * @param type      TypeID object to be used by default
         */
        virtual BasicModelForPCE& setDefaultObs( const SatID::SatelliteSystem& sys
                                                            = SatID::systemGPS,
                                                 const TypeID& type = TypeID::PC )
        {
            if(sys == SatID::systemGPS)
                defaultObsOfGPS = type;
            else if(sys == SatID::systemGalileo)
                defaultObsOfGAL = type;
            else if(sys == SatID::systemBDS)
                defaultObsOfBDS = type;

            return (*this);
        };


        /// Method to get use of TGD
        virtual bool getUseTGD(const SatID::SatelliteSystem sys = SatID::systemGPS)
        {
            if(sys == SatID::systemGPS)
                return useTGDOfGPS;
            else if(sys == SatID::systemGalileo)
                return useTGDOfGAL;
            else if(sys == SatID::systemBDS)
                return useTGDOfBDS;
        };

        /// Method to set use of TGD
        virtual BasicModelForPCE& setUseTGD(bool tgd,
                    const SatID::SatelliteSystem& sys = SatID::systemGPS)
        {
            if(sys == SatID::systemGPS)
                useTGDOfGPS = tgd;
            else if(sys == SatID::systemGalileo)
                useTGDOfGAL = tgd;
            else if(sys == SatID::systemBDS)
                useTGDOfBDS = tgd;

            return (*this);
        };


        /// Method to get Position to be used with GNSS data structures.
        virtual Position getNominalPosition() const
        { return nominalPos; };


        /** Method to set the station position to be used with GNSS
         *  data structures.
         *
         * @param pos     Position object to be used
         */
        virtual BasicModelForPCE& setNominalPosition(const Position& pos)
        { nominalPos = pos; return (*this); };


        /// Method to get a pointer to the default XvtStore<SatID> to be used
        /// with GNSS data structures.
        virtual XvtStore<SatID>* getSP3Store() const
        { return pSP3Store; };


        /** Method to set the default XvtStore<SatID> to be used with GNSS
         *  data structures.
         *
         * @param ephStore     XvtStore<SatID> object to be used by default
         */
        virtual BasicModelForPCE& setSP3Store(XvtStore<SatID>& ephStore)
        { pSP3Store = &ephStore; return (*this); };


        /// Method to get a pointer to the default XvtStore<SatID> to be used
        /// with GNSS data structures.
        virtual XvtStore<SatID>* getBCEStore() const
        { return pBCEStore; };


        /** Method to set the default XvtStore<SatID> to be used with GNSS
         *  data structures.
         *
         * @param ephStore     XvtStore<SatID> object to be used by default
         */
        virtual BasicModelForPCE& setBCEStore(XvtStore<SatID>& ephStore)
        { pBCEStore = &ephStore; return (*this); };


        /// Return a pointer to the MSCStore object currently in use.
        virtual MSCStore *getMSCStore(void) const
        { return pMSCStore; };


        /** Sets MSCStore object to be used.
         *
         * @param msc  MSCStore object.
         */
        virtual BasicModelForPCE& setMSCStore(MSCStore& msc)
        { pMSCStore = &msc; return (*this); };


        /// Get satellite clock.
        virtual satValueMap getSatClockSP3() const
        { return satClockSP3; };


        /// Get satellite clock.
        virtual satValueMap getSatClockBCE() const
        { return satClockBCE; };


        /// Return a string identifying this object.
        virtual std::string getClassName(void) const;


        /// Destructor.
        virtual ~BasicModelForPCE() {};


    protected:


        /// The elevation cut-off angle for accepted satellites.
        /// By default it is set to 10 degrees.
        double minElev;

        /// Station position
        Position nominalPos;

        /// Pointer to XvtStore<SatID> object
        XvtStore<SatID>* pSP3Store;

        /// Pointer to XvtStore<SatID> object
        XvtStore<SatID>* pBCEStore;

        /// Satellite clock
        satValueMap satClockSP3;
        satValueMap satClockBCE;

        /// Pointer to MSCStore object
        MSCStore* pMSCStore;


        /// Default observable to be used
        TypeID defaultObsOfGPS;
        TypeID defaultObsOfGAL;
        TypeID defaultObsOfBDS;


        /// Whether the TGD effect will be applied to C1 observable or not.
        bool useTGDOfGPS;
        bool useTGDOfGAL;
        bool useTGDOfBDS;


        /** Method to set the initial (a priori) position of receiver.
         * @return
         *  0 if OK
         *  -1 if problems arose
         */
        virtual int setInitialStaPosition( const double& aRx,
                                           const double& bRx,
                                           const double& cRx,
                           Position::CoordinateSystem s = Position::Cartesian,
                                         EllipsoidModel *ell = NULL,
                           ReferenceFrame frame = ReferenceFrame::Unknown );


        /// Method to set the initial (a priori) position of receiver.
        virtual int setInitialStaPosition(const Position& StaCoordinates);


        /// Method to set the initial (a priori) position of receiver.
        virtual int setInitialStaPosition();


        /// Method to get TGD corrections.
        virtual double getTGDCorrections( CommonTime Tr,
                                          const XvtStore<SatID>& Eph,
                                          SatID sat,
                                          TypeID type );


    }; // End of class 'BasicModelForPCE'

    //@}

}  // End of namespace gpstk

#endif   // BASICMODEL_FOR_PCE_HPP
