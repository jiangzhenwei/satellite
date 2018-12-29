#pragma ident "$Id$"

/**
 * @file BasicModel.cpp
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


#include "BasicModel.hpp"
#include "Rinex3EphemerisStore2.hpp"
#include "YDSTime.hpp"
#include "constants.hpp"

using namespace std;

namespace gpstk
{


      // Return a string identifying this object.
    std::string BasicModel::getClassName() const
    { return "BasicModel"; }



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
    BasicModel::BasicModel( const double& aSta,
                            const double& bSta,
                            const double& cSta,
                            Position::CoordinateSystem s,
                            EllipsoidModel *ell,
                            ReferenceFrame frame )
    {

        minElev = 10.0;
        nominalPos = Position(aSta, bSta, cSta, s, ell, frame);
        pEphStore = NULL;
        pMSCStore = NULL;
        defaultObsOfGPS = TypeID::C1;
        defaultObsOfGAL = TypeID::C1;
        defaultObsOfBDS = TypeID::C2;
        useTGDOfGPS = false;
        useTGDOfGAL = false;
        useTGDOfBDS = false;

    }  // End of 'BasicModel::BasicModel()'


      // Explicit constructor, taking as input a Position object
      // containing reference station coordinates.
    BasicModel::BasicModel(const Position& staPos)
    {

        minElev = 10.0;
        nominalPos = staPos;
        pEphStore = NULL;
        pMSCStore = NULL;
        defaultObsOfGPS = TypeID::C1;
        defaultObsOfGAL = TypeID::C1;
        defaultObsOfBDS = TypeID::C2;
        useTGDOfGPS = false;
        useTGDOfGAL = false;
        useTGDOfBDS = false;

    }  // End of 'BasicModel::BasicModel()'



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
    BasicModel::BasicModel( const Position& StaPos,
                            XvtStore<SatID>& ephStore,
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
        pEphStore = &ephStore;
        pMSCStore = &mscStore;
        defaultObsOfGPS = dObsOfGPS;
        defaultObsOfGAL = dObsOfGAL;
        defaultObsOfBDS = dObsOfBDS;
        useTGDOfGPS = applyTGDOfGPS;
        useTGDOfGAL = applyTGDOfGAL;
        useTGDOfBDS = applyTGDOfBDS;

    }  // End of 'BasicModel::BasicModel()'



      /* Return a satTypeValueMap object, adding the new data generated when
       * calling a modeling object.
       *
       * @param time      Epoch.
       * @param gData     Data object holding the data.
       */
    satTypeValueMap& BasicModel::Process( const CommonTime& time,
                                          satTypeValueMap& gData )
        throw(ProcessingException)
    {

        try
        {
            SatIDSet satRejectedSet;

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


                CommonTime tt;
                CommonTime transmit( time );

                Xvt svPosVel;

                double tof(0.0);
                double wt(0.0);
                double rawrange(0.0);
                double relativity(0.0);
                double svclkbias(0.0), svclkdrift(0.0);
                Triple cosines;
                double elevation(0.0), azimuth(0.0);
                double elevationGeodetic(0.0), azimuthGeodetic(0.0);

                try
                {
                    transmit -= obs/C_MPS;
                    tt = transmit;

                    for(int i=0; i<2; i++)
                    {
                        svPosVel = pEphStore->getXvt(sat,tt);

                        tt = transmit;

                        tt -= (svPosVel.clkbias + svPosVel.relcorr);
                    }

                    // rotate from ECF at transmit time to ECF at receive time
                    tof = RSS(svPosVel.x[0] - nominalPos.X(),
                              svPosVel.x[1] - nominalPos.Y(),
                              svPosVel.x[2] - nominalPos.Z())/C_MPS;
                    wt = OMEGA_EARTH * tof;

                    double sx(0.0), sy(0.0);

                    sx = +std::cos(wt)*svPosVel.x[0] + std::sin(wt)*svPosVel.x[1];
                    sy = -std::sin(wt)*svPosVel.x[0] + std::cos(wt)*svPosVel.x[1];

                    svPosVel.x[0] = sx;
                    svPosVel.x[1] = sy;

                    sx = +std::cos(wt)*svPosVel.v[0] + std::sin(wt)*svPosVel.v[1];
                    sy = -std::sin(wt)*svPosVel.v[0] + std::cos(wt)*svPosVel.v[1];

                    svPosVel.v[0] = sx;
                    svPosVel.v[1] = sy;

                    // rawrange
                    rawrange = RSS(svPosVel.x[0] - nominalPos.X(),
                                   svPosVel.x[1] - nominalPos.Y(),
                                   svPosVel.x[2] - nominalPos.Z());

                    // relativity
                    relativity = svPosVel.computeRelativityCorrection()*C_MPS;

                    // clock bias, clock drift
                    svclkbias = svPosVel.clkbias*C_MPS;
                    svclkdrift = svPosVel.clkdrift*C_MPS;

                    // partials
                    cosines[0] = (nominalPos.X() - svPosVel.x[0]) / rawrange;
                    cosines[1] = (nominalPos.Y() - svPosVel.x[1]) / rawrange;
                    cosines[2] = (nominalPos.Z() - svPosVel.x[2]) / rawrange;

                    // elevation, azimuth
                    Position SV( svPosVel );
                    elevation = nominalPos.elevation(SV);
                    azimuth = nominalPos.azimuth(SV);
                    elevationGeodetic = nominalPos.elevationGeodetic(SV);
                    azimuthGeodetic = nominalPos.azimuthGeodetic(SV);

                }
                catch(InvalidRequest& e)
                {
                    // If some problem appears, then schedule this satellite
                    // for removal
                    satRejectedSet.insert( sat );
                    continue;    // Skip this SV if problems arise
                }

                // Let's test if satellite has enough elevation over horizon
                if ( nominalPos.elevationGeodetic(svPosVel) < minElev )
                {
                    // Mark this satellite if it doesn't have enough elevation
                    satRejectedSet.insert( sat );
                    continue;
                }

                if(satClockMap.find(sat) == satClockMap.end())
                {
                    satClockMap[sat] = svclkbias;
                }

                // rho, relativity, elevation, azimuth
                (*it).second[TypeID::rho] = rawrange;
                (*it).second[TypeID::relativity] = -relativity;
                (*it).second[TypeID::elevation] = elevationGeodetic;
                (*it).second[TypeID::azimuth] = azimuthGeodetic;

                // Let's insert satellite position at transmit time
                (*it).second[TypeID::satXECF] = svPosVel.x[0];
                (*it).second[TypeID::satYECF] = svPosVel.x[1];
                (*it).second[TypeID::satZECF] = svPosVel.x[2];

                // Let's insert satellite velocity at transmit time
                (*it).second[TypeID::satVXECF] = svPosVel.v[0];
                (*it).second[TypeID::satVYECF] = svPosVel.v[1];
                (*it).second[TypeID::satVZECF] = svPosVel.v[2];

                // Let's insert satellite clock bias at transmit time
                (*it).second[TypeID::cdtSat] = svclkbias;
                (*it).second[TypeID::cdtSatDot] = svclkdrift;

                // Let's insert partials for station clock bias at receive time
                (*it).second[TypeID::dcdtSta] = 1.0;

                // Let's insert station position at receive time
                (*it).second[TypeID::staXECF] = nominalPos.X();
                (*it).second[TypeID::staYECF] = nominalPos.Y();
                (*it).second[TypeID::staZECF] = nominalPos.Z();

                // Let's insert partials for satellite position at transmit time
                (*it).second[TypeID::dSatX] = -cosines[0];
                (*it).second[TypeID::dSatY] = -cosines[1];
                (*it).second[TypeID::dSatZ] = -cosines[2];

                // Let's insert partials for station position at receive time
                (*it).second[TypeID::dStaX] = cosines[0];
                (*it).second[TypeID::dStaY] = cosines[1];
                (*it).second[TypeID::dStaZ] = cosines[2];

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

    }  // End of method 'BasicModel::Process()'


     /** Return a gnssDataMap object, adding the new data generated when
      *  calling a modeling object.
      *
      * @param gData    Data object holding the data.
      */
    gnssDataMap& BasicModel::Process(gnssDataMap& gData)
        throw(ProcessingException)
    {

        SourceIDSet sourceRejectedSet;

        satClockMap.clear();

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
                string station( source.sourceName );

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

    }  // End of method 'BasicModel::Process()'

}  // End of namespace gpstk
