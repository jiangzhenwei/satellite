#pragma ident "$Id$"

/**
 * @file BasicModelForPCE.cpp
 * This is a class to compute the basic parts of a GNSS model, i.e.:
 * Geometric distance, relativity correction, satellite position and
 * velocity at transmission time, satellite elevation and azimuth, etc.
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
//  Dagoberto Salazar - gAGE ( http://www.gage.es ). 2007, 2008, 2009, 2011
//
//============================================================================


#include "BasicModelForPCE.hpp"
#include "RinexEphemerisStore.hpp"
#include "Rinex3EphemerisStore2.hpp"
#include "YDSTime.hpp"
#include "constants.hpp"

using namespace std;

namespace gpstk
{


      // Return a string identifying this object.
   std::string BasicModelForPCE::getClassName() const
   { return "BasicModelForPCE"; }



      /* Explicit constructor taking as input reference
       * station coordinates.
       *
       * Those coordinates may be Cartesian (X, Y, Z in meters) or Geodetic
       * (Latitude, Longitude, Altitude), but defaults to Cartesian.
       *
       * Also, a pointer to GeoidModel may be specified, but default is
       * NULL (in which case WGS84 values will be used).
       *
       * @param aRx   first coordinate [ X(m), or latitude (degrees N) ]
       * @param bRx   second coordinate [ Y(m), or longitude (degrees E) ]
       * @param cRx   third coordinate [ Z, height above ellipsoid or
       *              radius, in meters ]
       * @param s     coordinate system (default is Cartesian, may be set
       *              to Geodetic).
       * @param ell   pointer to EllipsoidModel
       * @param frame Reference frame associated with this position
       */
    BasicModelForPCE::BasicModelForPCE( const double& aSta,
                              const double& bSta,
                              const double& cSta,
                              Position::CoordinateSystem s,
                              EllipsoidModel *ell,
                              ReferenceFrame frame )
    {

        minElev = 10.0;
        nominalPos = Position(aSta, bSta, cSta, s, ell, frame);
        pSP3Store = NULL;
        pBCEStore = NULL;
        pMSCStore = NULL;
        defaultObsOfGPS = TypeID::C1;
        defaultObsOfGAL = TypeID::C1;
        defaultObsOfBDS = TypeID::C2;
        useTGDOfGPS = false;
        useTGDOfGAL = false;
        useTGDOfBDS = false;

    }  // End of 'BasicModelForPCE::BasicModelForPCE()'


      // Explicit constructor, taking as input a Position object
      // containing reference station coordinates.
    BasicModelForPCE::BasicModelForPCE(const Position& staPos)
    {

        minElev = 10.0;
        nominalPos = staPos;
        pSP3Store = NULL;
        pBCEStore = NULL;
        pMSCStore = NULL;
        defaultObsOfGPS = TypeID::C1;
        defaultObsOfGAL = TypeID::C1;
        defaultObsOfBDS = TypeID::C2;
        useTGDOfGPS = false;
        useTGDOfGAL = false;
        useTGDOfBDS = false;

    }  // End of 'BasicModelForPCE::BasicModelForPCE()'



      /* Explicit constructor, taking as input reference station
       * coordinates, ephemeris to be used, default observable
       * and whether TGD will be computed or not.
       *
       * @param RxCoordinates Reference station coordinates.
       * @param dEphemeris    EphemerisStore object to be used by default.
       * @param dObservable   Observable type to be used by default.
       * @param applyTGD      Whether or not C1 observable will be
       *                      corrected from TGD effect or not.
       *
       */
    BasicModelForPCE::BasicModelForPCE( const Position& StaPos,
                              XvtStore<SatID>& sp3Store,
                              XvtStore<SatID>& bceStore,
                              MSCStore& mscStore,
                              const TypeID& dObsOfGPS,
                              const TypeID& dObsOfGAL,
                              const TypeID& dObsOfBDS,
                              const bool& applyTGDOfGPS,
                              const bool& applyTGDOfGAL,
                              const bool& applyTGDOfBDS )
    {

        minElev = 10.0;
        nominalPos = StaPos;
        pSP3Store = &sp3Store;
        pBCEStore = &bceStore;
        pMSCStore = &mscStore;
        defaultObsOfGPS = dObsOfGPS;
        defaultObsOfGAL = dObsOfGAL;
        defaultObsOfBDS = dObsOfBDS;
        useTGDOfGPS = applyTGDOfGPS;
        useTGDOfGAL = applyTGDOfGAL;
        useTGDOfBDS = applyTGDOfBDS;

    }  // End of 'BasicModelForPCE::BasicModelForPCE()'



      /* Return a satTypeValueMap object, adding the new data generated when
       * calling a modeling object.
       *
       * @param time      Epoch.
       * @param gData     Data object holding the data.
       */
    satTypeValueMap& BasicModelForPCE::Process( const CommonTime& time,
                                           satTypeValueMap& gData )
        throw(ProcessingException)
    {

        try
        {
            // Satellites to be rejected
            SatIDSet satRejectedSet;

            // Type
            TypeID defaultObs;

            // Loop through all the satellites
            for(satTypeValueMap::iterator it = gData.begin();
                it != gData.end();
                ++it)
            {
                SatID sat( it->first );

                // Scalar to hold temporal value
                double obs(0.0);

                if(sat.system == SatID::systemGPS)
                {
                    defaultObs = defaultObsOfGPS;
                }
                else if(sat.system == SatID::systemGalileo)
                {
                    defaultObs = defaultObsOfGAL;
                }
                else if(sat.system == SatID::systemBDS)
                {
                    defaultObs = defaultObsOfBDS;
                }

                obs = (*it).second(defaultObs);

                // A lot of the work is done by a CorrectedEphemerisRange object
                CorrectedEphemerisRange cerangeSP3;
                CorrectedEphemerisRange cerangeBCE;

                try
                {
                    // Compute most of the parameters
                    cerangeSP3.ComputeAtTransmitTime( time,
                                                      obs,
                                                      nominalPos,
                                                      sat,
                                                      *pSP3Store );
                    cerangeBCE.ComputeAtTransmitTime( time,
                                                      obs,
                                                      nominalPos,
                                                      sat,
                                                      *pBCEStore );
                }
                catch(InvalidRequest& e)
                {
                    // If some problem appears, then schedule this satellite
                    // for removal
                    satRejectedSet.insert( sat );
                    continue;    // Skip this SV if problems arise
                }

                // Let's test if satellite has enough elevation over horizon
                if ( nominalPos.elevationGeodetic(cerangeSP3.svPosVel) < minElev )
                {
                    // Mark this satellite if it doesn't have enough elevation
                    satRejectedSet.insert( sat );
                    continue;
                }


                if(satClockSP3.find(sat) == satClockSP3.end())
                {
                    double clock(cerangeSP3.svclkbias);
                    satClockSP3[sat] = clock;
                }

                if(satClockBCE.find(sat) == satClockBCE.end())
                {
                    double clock(cerangeBCE.svclkbias);
                    satClockBCE[sat] = clock;
                }

                // Computing Total Group Delay (TGD - meters), if possible
                double tempTGD(0.0);

                // coefficients of satellite position in ITRS
                (*it).second[TypeID::dSatX] = -cerangeSP3.cosines[0];
                (*it).second[TypeID::dSatY] = -cerangeSP3.cosines[1];
                (*it).second[TypeID::dSatZ] = -cerangeSP3.cosines[2];

                // coefficients of station position in ITRS
                (*it).second[TypeID::dStaX] = cerangeSP3.cosines[0];
                (*it).second[TypeID::dStaY] = cerangeSP3.cosines[1];
                (*it).second[TypeID::dStaZ] = cerangeSP3.cosines[2];

                // Now we have to add the new values to the data structure
                (*it).second[TypeID::rho] = cerangeSP3.rawrange;
                (*it).second[TypeID::relativity] = -cerangeSP3.relativity;
                (*it).second[TypeID::elevation] = cerangeSP3.elevationGeodetic;
                (*it).second[TypeID::azimuth] = cerangeSP3.azimuthGeodetic;

                // Let's insert satellite position at transmission time
                (*it).second[TypeID::satXECF] = cerangeSP3.svPosVel.x[0];
                (*it).second[TypeID::satYECF] = cerangeSP3.svPosVel.x[1];
                (*it).second[TypeID::satZECF] = cerangeSP3.svPosVel.x[2];

                // Let's insert satellite velocity at transmission time
                (*it).second[TypeID::satVXECF] = cerangeSP3.svPosVel.v[0];
                (*it).second[TypeID::satVYECF] = cerangeSP3.svPosVel.v[1];
                (*it).second[TypeID::satVZECF] = cerangeSP3.svPosVel.v[2];

                // Let's insert satellite clock bias at transmission time
                (*it).second[TypeID::cdtSat] = cerangeBCE.svclkbias;
                (*it).second[TypeID::cdtSatDot] = cerangeBCE.svclkdrift;

                // Let's insert station position
                (*it).second[TypeID::staXECF] = nominalPos.X();
                (*it).second[TypeID::staYECF] = nominalPos.Y();
                (*it).second[TypeID::staZECF] = nominalPos.Z();

            } // End of loop for(it = gData.begin()...

            // Remove satellites with missing data
            gData.removeSatID(satRejectedSet);

            return gData;

        }   // End of try...
        catch(Exception& u)
        {
            // Throw an exception if something unexpected happens
            ProcessingException e( getClassName() + ":" + u.what() );
            GPSTK_THROW(e);
        }

    }  // End of method 'BasicModelForPCE::Process()'


     /** Return a gnssDataMap object, adding the new data generated when
      *  calling a modeling object.
      *
      * @param gData    Data object holding the data.
      */
    gnssDataMap& BasicModelForPCE::Process(gnssDataMap& gData)
        throw(ProcessingException)
    {
        SourceIDSet sourceRejectedSet;

        satClockSP3.clear();
        satClockBCE.clear();

        for(gnssDataMap::iterator gdmIt = gData.begin();
            gdmIt != gData.end();
            ++gdmIt)
        {
            CommonTime epoch( gdmIt->first );
            epoch.setTimeSystem( TimeSystem::Unknown );

            for(sourceDataMap::iterator sdmIt = gdmIt->second.begin();
                sdmIt != gdmIt->second.end();
                ++sdmIt)
            {
                SourceID source( sdmIt->first );
                string station( source.sourceName.substr(0,4) );

                if(pMSCStore == NULL)
                {
                    sourceRejectedSet.insert( source );
                    continue;
                }

                MSCData mscData;

                try
                {
                    mscData = pMSCStore->findMSC(station,epoch);
                }
                catch(...)
                {
                    sourceRejectedSet.insert( source );
                    continue;
                }

                nominalPos = mscData.coordinates;

                Process( gdmIt->first, sdmIt->second );
            }
        }

        gData.removeSourceID( sourceRejectedSet );

        return gData;

    }  // End of method 'BasicModelForPCE::Process()'



      /* Method to set the initial (a priori) position of receiver.
       * @return
       *  0 if OK
       *  -1 if problems arose
       */
    int BasicModelForPCE::setInitialStaPosition( const double& aSta,
                                                 const double& bSta,
                                                 const double& cSta,
                                                 Position::CoordinateSystem s,
                                                 EllipsoidModel *ell,
                                                 ReferenceFrame frame )
    {

        try
        {
            Position stapos( aSta, bSta, cSta, s, ell, frame );
            setInitialStaPosition(stapos);
            return 0;
        }
        catch(GeometryException& e)
        {
            return -1;
        }

    }  // End of method 'BasicModelForPCE::setInitialStaPosition()'



      // Method to set the initial (a priori) position of receiver.
    int BasicModelForPCE::setInitialStaPosition(const Position& StaCoordinates)
    {

        try
        {
            nominalPos = StaCoordinates;
            return 0;
        }
        catch(GeometryException& e)
        {
            return -1;
        }

    }  // End of method 'BasicModelForPCE::setInitialStaPosition()'



      // Method to set the initial (a priori) position of receiver.
    int BasicModelForPCE::setInitialStaPosition()
    {
        try
        {
            Position stapos(0.0, 0.0, 0.0, Position::Cartesian, NULL);
            setInitialStaPosition(stapos);
            return 0;
        }
        catch(GeometryException& e)
        {
            return -1;
        }

    }  // End of method 'BasicModelForPCE::setInitialStaPosition()'



      // Method to get TGD corrections.
    double BasicModelForPCE::getTGDCorrections( CommonTime Tr,
                                                const XvtStore<SatID>& Eph,
                                                SatID sat,
                                                TypeID type )
    {

        try
        {
            const Rinex3EphemerisStore2& bce =
                            dynamic_cast<const Rinex3EphemerisStore2&>(Eph);

            if(sat.system == SatID::systemGPS)
            {
                const GPSEphemerisStore& gps = bce.getGPSEphemerisStore();

                Tr.setTimeSystem(TimeSystem::GPS);
                double tgd = gps.findEphemeris(sat,Tr).Tgd * C_MPS;

                if(type == TypeID::C1)        return ( 1.0 * tgd );
                else if(type == TypeID::C2)   return ( GAMMA_GPS_L1L2 * tgd );
                else if(type == TypeID::PC)   return ( 0.0 );
            }
            else if(sat.system == SatID::systemGalileo)
            {
                const GalEphemerisStore& gal = bce.getGalEphemerisStore();

                Tr.setTimeSystem(TimeSystem::GAL);
                double tgda = gal.findEphemeris(sat,Tr).Tgda * C_MPS;

                if(type == TypeID::C1)      return ( 1.0 * tgda );
                else if(type == TypeID::C5) return ( GAMMA_GAL_L1L5 * tgda );
                else if(type == TypeID::PC) return ( 0.0 );
            }
            else if(sat.system == SatID::systemBDS)
            {
                const BDSEphemerisStore& bds = bce.getBDSEphemerisStore();

                Tr.setTimeSystem(TimeSystem::BDT);
                double tgd13 = bds.findEphemeris(sat,Tr).Tgd13 * C_MPS;
                double tgd23 = bds.findEphemeris(sat,Tr).Tgd23 * C_MPS;

                if(type == TypeID::C2)      return ( 1.0 * tgd13 );
                else if(type == TypeID::C7) return ( 1.0 * tgd23 );
                else if(type == TypeID::PC) return ( 0.0 );
            }
        }
        catch(...)
        {
            cerr << "Exception " << sat << ' ' << type << endl;
            return 0.0;
        }

    }  // End of method 'BasicModelForPCE::getTGDCorrections()'

}  // End of namespace gpstk
