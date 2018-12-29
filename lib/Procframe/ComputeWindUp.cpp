#pragma ident "$Id$"

/**
 * @file ComputeWindUp.cpp
 * This class computes the wind-up effect on the phase observables, in radians.
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


#include "ComputeWindUp.hpp"

using namespace std;

namespace gpstk
{

    // Return a string identifying this object.
    std::string ComputeWindUp::getClassName() const
    { return "ComputeWindUp"; }



    /** Return a satTypeValueMap object, adding the new data generated when
     *  calling this object.
     *
     * @param time      Epoch corresponding to the data.
     * @param gData     Data object holding the data.
     */
    satTypeValueMap& ComputeWindUp::Process( const CommonTime& time,
                                              satTypeValueMap& gData )
        throw(ProcessingException)
    {
        try
        {
            if(pRefSys == NULL)
            {
                ProcessingException e("Pointer to ReferenceSystem NULL!");
                GPSTK_THROW(e);
            }

            if(pSolSys == NULL)
            {
                ProcessingException e("Pointer to SolarSystem NULL!");
                GPSTK_THROW(e);
            }

            if(pEphStore == NULL)
            {
                ProcessingException e("Pointer to EphemerisStore NULL!");
                GPSTK_THROW(e);
            }

            CommonTime utc( pRefSys->GPS2UTC(time) );

            Matrix<double> t2cRaw( pRefSys->T2CMatrix(utc) );
            Matrix<double> t2cDot( pRefSys->dT2CMatrix(utc) );

            Matrix<double> c2tRaw( transpose( t2cRaw ) );

            CommonTime tt( pRefSys->GPS2TT(time) );
            double jd_tt( JulianDate(tt).jd );

            SolarSystem::Planet center(SolarSystem::Earth);
            SolarSystem::Planet target(SolarSystem::Sun);

            // Sun pos and vel in ECI, unit: km, km/day
            double rv_sun[6] = {0.0};
            pSolSys->computeState(jd_tt, target, center, rv_sun);

            // Sun pos in ECI, unit: m
            Vector<double> rSunECI(3,0.0);
            rSunECI(0) = rv_sun[0];
            rSunECI(1) = rv_sun[1];
            rSunECI(2) = rv_sun[2];
            rSunECI *= 1.0e+3;

            Vector<double> rSunECF( t2cRaw*rSunECI );

            // vector from Earth to Station
            Vector<double> rStaECF(3,0.0);
            rStaECF(0) = nominalPos.X();
            rStaECF(1) = nominalPos.Y();
            rStaECF(2) = nominalPos.Z();

            // Sat pos in ECF, unit: m
            Vector<double> rSatECF(3, 0.0);

            // Sat vel in ECF, unit: m/s
            Vector<double> vSatECF(3, 0.0);

            // Sat pos in ECI, unit: m
            Vector<double> rSatECI(3, 0.0);

            // Sat vel in ECI, unit: m/s
            Vector<double> vSatECI(3, 0.0);


            SatIDSet satRejectedSet;

            // Loop through all the satellites
            for( satTypeValueMap::iterator it = gData.begin();
                 it != gData.end();
                 ++it )
            {
                SatID sat( it->first );

                // First check if this satellite has previous arc information
                if( satPhaseData.find( sat ) == satPhaseData.end() )
                {
                    // If it doesn't have an entry, insert one
                    satPhaseData[ sat ].arcNum = 0.0;
                }

                // Then, check both if there is arc information, and if current
                // arc number is different from arc number in storage (which
                // means a cycle slip happened)
                if( it->second.find(TypeID::satArc) != it->second.end() &&
                    it->second(TypeID::satArc) != satPhaseData[ sat ].arcNum )
                {
                    // If different, update satellite arc in storage
                    satPhaseData[ sat ].arcNum = it->second(TypeID::satArc);

                    // Reset phase information
                    satPhaseData[ sat ].satPreviousPhase = 0.0;
                    satPhaseData[ sat ].staPreviousPhase = 0.0;
                }


                // Use ephemeris if satellite position is not already computed
                if( ( it->second.find(TypeID::satXECF) == it->second.end() ) ||
                    ( it->second.find(TypeID::satYECF) == it->second.end() ) ||
                    ( it->second.find(TypeID::satZECF) == it->second.end() ) ||
                    ( it->second.find(TypeID::satVXECF) == it->second.end() ) ||
                    ( it->second.find(TypeID::satVYECF) == it->second.end() ) ||
                    ( it->second.find(TypeID::satVZECF) == it->second.end() ) )
                {
                    try
                    {
                        rSatECF = pEphStore->getXvt(sat,time).x.toVector();
                        vSatECF = pEphStore->getXvt(sat,time).v.toVector();
                    }
                    catch(...)
                    {
                        satRejectedSet.insert( sat );
                        continue;
                    }
                }
                else
                {
                    // Get satellite position out of GDS
                    rSatECF[0] = it->second[TypeID::satXECF];
                    rSatECF[1] = it->second[TypeID::satYECF];
                    rSatECF[2] = it->second[TypeID::satZECF];

                    vSatECF[0] = it->second[TypeID::satVXECF];
                    vSatECF[1] = it->second[TypeID::satVYECF];
                    vSatECF[2] = it->second[TypeID::satVZECF];
                }

                rSatECI = t2cRaw * rSatECF;
//                vSatECI = t2cRaw * vSatECF + t2cDot * rSatECF;

                vSatECI(0) = vSatECF(0) - OMEGA_EARTH * rSatECF(1);
                vSatECI(1) = vSatECF(1) + OMEGA_EARTH * rSatECF(0);
                vSatECI(2) = vSatECF(2);


                // Get satellite rotation angle

                // unit vector of Sat to Earth
                Vector<double> rk( (-1.0)*normalize(rSatECF) );

                // unit vector of tangetial direction
                Vector<double> rt( normalize(vSatECI) );

                // modeled yaw attitude
                double yawDeg( satYawData[sat].modeled );
                double yawRad( yawDeg * DEG_TO_RAD );

                // quaternion rotate yaw with respect to unit_r
                Vector<double> qua(4,0.0);
                qua(0) = std::cos( yawRad/2.0 );
                qua(1) = rk(0) * std::sin( yawRad/2.0 );
                qua(2) = rk(1) * std::sin( yawRad/2.0 );
                qua(3) = rk(2) * std::sin( yawRad/2.0 );

                // matrix rotate yaw with respect to unit_r
                Matrix<double> mat(3,3,0.0);
                mat(0,0) = qua(0)*qua(0) + qua(1)*qua(1)
                         - qua(2)*qua(2) - qua(3)*qua(3);
                mat(0,1) = 2.0 * (qua(1)*qua(2) - qua(0)*qua(3));
                mat(0,2) = 2.0 * (qua(1)*qua(3) + qua(0)*qua(2));
                mat(1,0) = 2.0 * (qua(1)*qua(2) + qua(0)*qua(3));
                mat(1,1) = qua(0)*qua(0) - qua(1)*qua(1)
                         + qua(2)*qua(2) - qua(3)*qua(3);
                mat(1,2) = 2.0 * (qua(2)*qua(3) - qua(0)*qua(1));
                mat(2,0) = 2.0 * (qua(1)*qua(3) - qua(0)*qua(2));
                mat(2,1) = 2.0 * (qua(2)*qua(3) + qua(0)*qua(1));
                mat(2,2) = qua(0)*qua(0) - qua(1)*qua(1)
                         - qua(2)*qua(2) + qua(3)*qua(3);

                Vector<double> ri( normalize(mat*rt) );

                Vector<double> rj( normalize( cross(rk,ri) ) );


                // unit vector from Station to Sat
                Vector<double> rrho( normalize(rSatECF-rStaECF) );

                // Projection of "rk" vector to line of sight vector (rrho)
                double zk( dot(rrho,rk) );

                // Get a vector without components on rk (i.e., belonging
                // to ri, rj plane)
                Vector<double> dpp(rrho - zk*rk);

                // Compute dpp components in ri, rj plane
                double xk( dot(dpp,ri) );
                double yk( dot(dpp,rj) );

                // Compute satellite rotation angle, in radians
                double alpha1( std::atan2(yk,xk) );


                // Get station rotation angle

                // unit vector from Station to Earth
                rk = (-1.0)*normalize(rStaECF);

                // Let's define a NORTH unit vector in the Up, East, North
                // (UEN) topocentric reference frame
                Vector<double> delta(3,0.0);
                delta(2) = 1.0;

                double lat( nominalPos.geodeticLatitude() );
                double lon( nominalPos.longitude() );

                double slat( std::sin(lat*DEG_TO_RAD) );
                double clat( std::cos(lat*DEG_TO_RAD) );
                double slon( std::sin(lon*DEG_TO_RAD) );
                double clon( std::cos(lon*DEG_TO_RAD) );

                Matrix<double> rlat(3,3,0.0);
                rlat(0,0) = +clat; rlat(0,2) = -slat;
                rlat(1,1) = 1.0;
                rlat(2,0) = +slat; rlat(2,2) = +clat;

                Matrix<double> rlon(3,3,0.0);
                rlon(0,0) = +clon; rlon(0,1) = -slon;
                rlon(1,0) = +slon; rlat(1,1) = +clon;
                rlon(2,2) = 1.0;

                // Rotate delta to XYZ reference frame
                delta = rlon*(rlat*delta);


                // Computation of reference frame unitary vectors for station
                // rj = rk x delta
                rj = normalize( cross(rk,delta) );

                // ri = rj x rk
                ri = normalize( cross(rj,rk) );

                // Projection of "rk" vector to line of sight vector (rrho)
                zk = dot(rrho,rk);

                // Get a vector without components on rk (i.e., belonging
                // to ri, rj plane)
                dpp = rrho - zk*rk;

                // Compute dpp components in ri, rj plane
                xk = dot(dpp,ri);
                yk = dot(dpp,rj);

                // Compute station rotation angle, in radians
                double alpha2( std::atan2(yk,xk) );

                // Compute wind up effect in radians
                double windUp(0.0);

                alpha1 = alpha1 + windUp;

                double da1(alpha1-satPhaseData[sat].satPreviousPhase);
                double sda1( std::sin(da1) ), cda1( std::cos(da1) );

                double da2(alpha2-satPhaseData[sat].staPreviousPhase);
                double sda2( std::sin(da2) ), cda2( std::cos(da2) );

                // Let's avoid problems when passing from 359 to 0 degrees.
                satPhaseData[sat].satPreviousPhase += std::atan2(sda1, cda1);

                satPhaseData[sat].staPreviousPhase += std::atan2(sda2, cda2);

                windUp = satPhaseData[sat].satPreviousPhase -
                         satPhaseData[sat].staPreviousPhase;

                // Let's get wind-up value in radians, and insert it
                // into GNSS data structure.
                it->second[TypeID::windUp] = windUp;

//                cout << setw(10) << windUp << endl;

            }  // End of 'for (it = gData.begin(); it != gData.end(); ++it)'

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

    }  // End of method 'ComputeWindUp::Process()'



    /** Return a gnssDataMap object, adding the new data generated
     *  when calling this object.
     *
     * @param gData    Data object holding the data.
     */
    gnssDataMap& ComputeWindUp::Process(gnssDataMap& gData)
        throw(ProcessingException)
    {
        SourceIDSet sourceRejectedSet;

        if(pMSCStore == NULL)
        {
            ProcessingException e("Pointer to ReferenceSystem NULL!");
            GPSTK_THROW(e);
        }

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

                SatPhaseDataMap::iterator spdmIt = satPhaseDataMap.find(source);

                if(spdmIt != satPhaseDataMap.end())
                {
                    satPhaseData = spdmIt->second;
                }
                else
                {
                    satPhaseData = SatPhaseData();
                }

                Process( gdmIt->first, sdmIt->second );

                satPhaseDataMap[source] = satPhaseData;
//                break;
            }
//            break;
        }

        gData.removeSourceID( sourceRejectedSet );

        return gData;

    }  // End of method 'ComputeWindUp::Process()'

}  // End of namespace gpstk
