#pragma ident "$Id$"

/**
 * @file ComputeSatPCenter.cpp
 * This class computes the satellite antenna phase correction, in meters.
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
//  Dagoberto Salazar - gAGE ( http://www.gage.es ). 2008, 2009, 2011
//
//============================================================================


#include "ComputeSatPCenter.hpp"

using namespace std;

namespace gpstk
{

    // Return a string identifying this object.
    std::string ComputeSatPCenter::getClassName() const
    { return "ComputeSatPCenter"; }


    /** Return a satTypeValueMap object, adding the new data generated when
     *  calling this object.
     *
     * @param time      Epoch corresponding to the data.
     * @param gData     Data object holding the data.
     */
    satTypeValueMap& ComputeSatPCenter::Process( const CommonTime& time,
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

            if(pAntexReader == NULL)
            {
                ProcessingException e("Pointer to AntexReader NULL!");
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

            // Sun pos in ECF, unit: m
            Vector<double> rSunECF( c2tRaw*rSunECI );

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
            for(satTypeValueMap::iterator it = gData.begin();
                it != gData.end();
                ++it)
            {
                SatID sat( it->first );

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


                // Let's get the satellite antenna phase correction value in
                // meters, and insert it in the GNSS data structure.

                double satPCenterL1(0.0);
                double satPCenterL2(0.0);
                double satPCenterL3(0.0);
                double satPCenterL5(0.0);
                double satPCenterL6(0.0);
                double satPCenterL7(0.0);
                double satPCenterL8(0.0);
                double satPCenterL9(0.0);


                // unit vector of Sat to Earth
                Vector<double> rk( (-1.0)*normalize(rSatECF) );

//                // unit vector of normal direction
//                Vector<double> rn( normalize( cross(rSatECI,vSatECI) ) );
//
//                // unit vector of tangetial direction
//                Vector<double> rt( normalize( cross(rn,rSatECI) ) );

                Vector<double> rt( normalize(vSatECI) );

                // modeled yaw attitude
                double yawDeg( satYawData[sat].modeled );
                double yawRad( yawDeg * DEG_TO_RAD );

                // quaternion, rotate yaw with respect to unit_r
                Vector<double> qua(4,0.0);
                qua(0) = std::cos( yawRad/2.0 );
                qua(1) = rk(0) * std::sin( yawRad/2.0 );
                qua(2) = rk(1) * std::sin( yawRad/2.0 );
                qua(3) = rk(2) * std::sin( yawRad/2.0 );

                // matrix, rotate yaw with respect to unit_r
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

//                cout << setw(10) << rj(0)
//                     << setw(10) << rj(1)
//                     << setw(10) << rj(2)
//                     << endl;
//
//                break;


                // unit vector from Station to Sat
                Vector<double> rrho( normalize(rSatECF-rStaECF) );


                // We will need the elevation, in degrees. It is found using
                // dot product and the corresponding unitary angles

                double temp( -dot(rrho,rk) );
                if(temp > +1.0) temp = +1.0;
                if(temp < -1.0) temp = -1.0;

                double nadir( std::acos(temp) * RAD_TO_DEG );

                double elev( 0.0 );

                Triple pco(0.0, 0.0, 0.0);
                Triple pcv(0.0, 0.0, 0.0);

                Vector<double> satAntenna(3,0.0);

                // Get satellite information
                if( sat.system == SatID::systemGPS )
                {
                    std::stringstream satstr;
                    satstr << "G";
                    if( sat.id < 10 ) satstr << "0";
                    satstr << sat.id;

                    // Get satellite antenna information out of AntexReader object
                    Antenna antenna( pAntexReader->getAntenna( satstr.str(), time ) );

                    double zen2( antenna.getZen2() );

                    nadir = (nadir>zen2) ? zen2 : nadir;

                    elev = 90.0 - nadir;

                    try
                    {
                        // Get antenna eccentricity for frequency "G01" (L1), in
                        // satellite reference system.
                        // NOTE: It is NOT in ECF, it is in UEN!!!
                        pco = antenna.getAntennaEccentricity(Antenna::G01);

                        // Now, get the phase center variation.
                        pcv = antenna.getAntennaPCVariation(Antenna::G01, elev);

                        // Change PCO to ECF
                        satAntenna = pco[2]*ri + pco[1]*rj + pco[0]*rk;

                        // Projection of "satAntenna" vector to line of sight vector rrho
                        // This correction is interpreted as an "advance" in the signal,
                        // instead of a delay. Therefore, it has negative sign
                        satPCenterL1 = dot(rrho,satAntenna) + pcv[0];

//                        cout << setw(10) << satPCenterL1 << endl;
                    }
                    catch(...)
                    {
                        satPCenterL1 = 0.0;
                    }

                    try
                    {
                        // Get antenna eccentricity for frequency "G02" (L2), in
                        // satellite reference system.
                        // NOTE: It is NOT in ECF, it is in UEN!!!
                        pco = antenna.getAntennaEccentricity(Antenna::G02);

                        // Now, get the phase center variation.
                        pcv = antenna.getAntennaPCVariation(Antenna::G02, elev);

                        // Change PCO to ECF
                        satAntenna = pco[2]*ri + pco[1]*rj + pco[0]*rk;

                        // Projection of "satAntenna" vector to line of sight vector rrho
                        // This correction is interpreted as an "advance" in the signal,
                        // instead of a delay. Therefore, it has negative sign
                        satPCenterL2 = dot(rrho,satAntenna) + pcv[0];
                    }
                    catch(...)
                    {
                        satPCenterL2 = 0.0;
                    }

                    try
                    {
                        // Get antenna eccentricity for frequency "G05" (L5), in
                        // satellite reference system.
                        // NOTE: It is NOT in ECF, it is in UEN!!!
                        pco = antenna.getAntennaEccentricity(Antenna::G02);

                        // Now, get the phase center variation.
                        pcv = antenna.getAntennaPCVariation(Antenna::G02, elev);

                        // Change PCO to ECF
                        satAntenna = pco[2]*ri + pco[1]*rj + pco[0]*rk;

                        // Projection of "satAntenna" vector to line of sight vector rrho
                        // This correction is interpreted as an "advance" in the signal,
                        // instead of a delay. Therefore, it has negative sign
                        satPCenterL5 = dot(rrho,satAntenna) + pcv[0];
                    }
                    catch(...)
                    {
                        satPCenterL5 = 0.0;
                    }
                }
                // Check if this satellite belongs to Galileo system
                else if( sat.system == SatID::systemGalileo )
                {
                    std::stringstream satstr;
                    satstr << "E";
                    if( sat.id < 10 ) satstr << "0";
                    satstr << sat.id;

                    // Get satellite antenna information out of AntexReader object
                    Antenna antenna( pAntexReader->getAntenna( satstr.str(), time ) );

                    double zen2( antenna.getZen2() );

                    nadir = (nadir>zen2) ? zen2 : nadir;

                    elev = 90.0 - nadir;

                    try
                    {
                        // Get antenna offset for frequency "E01" (Galileo), in
                        // satellite reference system.
                        // NOTE: It is NOT in ECF, it is in UEN!!!
                        pco = antenna.getAntennaEccentricity(Antenna::E01);

                        // Now, get the phase center variation.
                        pcv = antenna.getAntennaPCVariation(Antenna::E01, elev);

                        // Change PCO to ECF
                        satAntenna = pco[2]*ri + pco[1]*rj + pco[0]*rk;

                        // Project "satAntenna" vector to line of sight vector rrho
                        // This correction is interpreted as an "advance" in the signal,
                        // instead of a delay. Therefore, it has negative sign
                        satPCenterL1 = dot(rrho,satAntenna) + pcv[0];
                    }
                    catch(...)
                    {
                        satPCenterL1 = 0.0;
                    }

                    try
                    {
                        // Get antenna offset for frequency "E05" (Galileo), in
                        // satellite reference system.
                        // NOTE: It is NOT in ECF, it is in UEN!!!
                        pco = antenna.getAntennaEccentricity(Antenna::E05);

                        // Now, get the phase center variation.
                        pcv = antenna.getAntennaPCVariation(Antenna::E05, elev);

                        // Change PCO to ECF
                        satAntenna = pco[2]*ri + pco[1]*rj + pco[0]*rk;

                        // Project "satAntenna" vector to line of sight vector rrho
                        // This correction is interpreted as an "advance" in the signal,
                        // instead of a delay. Therefore, it has negative sign
                        satPCenterL5 = dot(rrho,satAntenna) + pcv[0];
                    }
                    catch(...)
                    {
                        satPCenterL5 = 0.0;
                    }

                    try
                    {
                        // Get antenna offset for frequency "E06" (Galileo), in
                        // satellite reference system.
                        // NOTE: It is NOT in ECF, it is in UEN!!!
                        pco = antenna.getAntennaEccentricity(Antenna::E06);

                        // Now, get the phase center variation.
                        pcv = antenna.getAntennaPCVariation(Antenna::E06, elev);

                        // Change PCO to ECF
                        satAntenna = pco[2]*ri + pco[1]*rj + pco[0]*rk;

                        // Project "satAntenna" vector to line of sight vector rrho
                        // This correction is interpreted as an "advance" in the signal,
                        // instead of a delay. Therefore, it has negative sign
                        satPCenterL6 = dot(rrho,satAntenna) + pcv[0];
                    }
                    catch(...)
                    {
                        satPCenterL6 = 0.0;
                    }

                    try
                    {
                        // Get antenna offset for frequency "E07" (Galileo), in
                        // satellite reference system.
                        // NOTE: It is NOT in ECF, it is in UEN!!!
                        pco = antenna.getAntennaEccentricity(Antenna::E07);

                        // Now, get the phase center variation.
                        pcv = antenna.getAntennaPCVariation(Antenna::E07, elev);

                        // Change PCO to ECF
                        satAntenna = pco[2]*ri + pco[1]*rj + pco[0]*rk;

                        // Project "satAntenna" vector to line of sight vector rrho
                        // This correction is interpreted as an "advance" in the signal,
                        // instead of a delay. Therefore, it has negative sign
                        satPCenterL7 = dot(rrho,satAntenna) + pcv[0];
                    }
                    catch(...)
                    {
                        satPCenterL7 = 0.0;
                    }

                    try
                    {
                        // Get antenna offset for frequency "E08" (Galileo), in
                        // satellite reference system.
                        // NOTE: It is NOT in ECF, it is in UEN!!!
                        pco = antenna.getAntennaEccentricity(Antenna::E08);

                        // Now, get the phase center variation.
                        pcv = antenna.getAntennaPCVariation(Antenna::E08, elev);

                        // Change PCO to ECF
                        satAntenna = pco[2]*ri + pco[1]*rj + pco[0]*rk;

                        // Project "satAntenna" vector to line of sight vector rrho
                        // This correction is interpreted as an "advance" in the signal,
                        // instead of a delay. Therefore, it has negative sign
                        satPCenterL8 = dot(rrho,satAntenna) + pcv[0];
                    }
                    catch(...)
                    {
                        satPCenterL8 = 0.0;
                    }
                }
                // Check if this satellite belongs to BeiDou system
                else if( sat.system == SatID::systemBDS )
                {
                    std::stringstream satstr;
                    satstr << "C";
                    if( sat.id < 10 ) satstr << "0";
                    satstr << sat.id;

                    if(sat.id > 30)
                    {
                        satPCenterL1 = 0.0;
                        satPCenterL2 = 0.0;
                        satPCenterL6 = 0.0;
                        satPCenterL7 = 0.0;
                        continue;
                    }

                    // Get satellite antenna information out of AntexReader object
                    Antenna antenna( pAntexReader->getAntenna( satstr.str(), time ) );

                    double zen2( antenna.getZen2() );

                    nadir = (nadir>zen2) ? zen2 : nadir;

                    elev = 90.0 - nadir;

                    try
                    {
                        // Get antenna offset for frequency "C01" (BDS), in
                        // satellite reference system.
                        // NOTE: It is NOT in ECF, it is in UEN!!!
                        pco = antenna.getAntennaEccentricity(Antenna::C01);

                        // Now, get the phase center variation.
                        pcv = antenna.getAntennaPCVariation(Antenna::C01, elev);

                        // Change PCO to ECF
                        satAntenna = pco[2]*ri + pco[1]*rj + pco[0]*rk;

                        // Project "satAntenna" vector to line of sight vector rrho
                        // This correction is interpreted as an "advance" in the signal,
                        // instead of a delay. Therefore, it has negative sign
                        satPCenterL1 = dot(rrho,satAntenna) + pcv[0];
                    }
                    catch(...)
                    {
                        satPCenterL1 = 0.0;
                    }

                    try
                    {
                        // Get antenna offset for frequency "C02" (BDS), in
                        // satellite reference system.
                        // NOTE: It is NOT in ECF, it is in UEN!!!
                        pco = antenna.getAntennaEccentricity(Antenna::C02);

                        // Now, get the phase center variation.
                        pcv = antenna.getAntennaPCVariation(Antenna::C02, elev);

                        // Change PCO to ECF
                        satAntenna = pco[2]*ri + pco[1]*rj + pco[0]*rk;

                        // Project "satAntenna" vector to line of sight vector rrho
                        // This correction is interpreted as an "advance" in the signal,
                        // instead of a delay. Therefore, it has negative sign
                        satPCenterL2 = dot(rrho,satAntenna) + pcv[0];
                    }
                    catch(...)
                    {
                        satPCenterL2 = 0.0;
                    }

                    try
                    {
                        // Get antenna offset for frequency "C06" (BDS), in
                        // satellite reference system.
                        // NOTE: It is NOT in ECF, it is in UEN!!!
                        pco = antenna.getAntennaEccentricity(Antenna::C06);

                        // Now, get the phase center variation.
                        pcv = antenna.getAntennaPCVariation(Antenna::C06, elev);

                        // Change PCO to ECF
                        satAntenna = pco[2]*ri + pco[1]*rj + pco[0]*rk;

                        // Project "satAntenna" vector to line of sight vector rrho
                        // This correction is interpreted as an "advance" in the signal,
                        // instead of a delay. Therefore, it has negative sign
                        satPCenterL6 = dot(rrho,satAntenna) + pcv[0];
                    }
                    catch(...)
                    {
                        satPCenterL6 = 0.0;
                    }

                    try
                    {
                        // Get antenna offset for frequency "C07" (BDS), in
                        // satellite reference system.
                        // NOTE: It is NOT in ECF, it is in UEN!!!
                        pco = antenna.getAntennaEccentricity(Antenna::C07);

                        // Now, get the phase center variation.
                        pcv = antenna.getAntennaPCVariation(Antenna::C07, elev);

                        // Change PCO to ECF
                        satAntenna = pco[2]*ri + pco[1]*rj + pco[0]*rk;

                        // Project "satAntenna" vector to line of sight vector rrho
                        // This correction is interpreted as an "advance" in the signal,
                        // instead of a delay. Therefore, it has negative sign
                        satPCenterL7 = dot(rrho,satAntenna) + pcv[0];
                    }
                    catch(...)
                    {
                        satPCenterL7 = 0.0;
                    }
                }


                // Find which observables are present, and then
                // apply corrections

                // Look for C1
                if(it->second.find(TypeID::C1) != it->second.end())
                {
                    it->second[TypeID::C1] -= satPCenterL1;
                }

                // Look for L1
                if( it->second.find(TypeID::L1) != it->second.end() )
                {
                    it->second[TypeID::L1] -= satPCenterL1;
                }

                // Look for C2
                if( it->second.find(TypeID::C2) != it->second.end() )
                {
                    it->second[TypeID::C2] -= satPCenterL2;
                }

                // Look for L2
                if( it->second.find(TypeID::L2) != it->second.end() )
                {
                    it->second[TypeID::L2] -= satPCenterL2;
                }

                // Look for C3
                if( it->second.find(TypeID::C3) != it->second.end() )
                {
                    it->second[TypeID::C3] -= satPCenterL3;
                }

                // Look for L3
                if( it->second.find(TypeID::L3) != it->second.end() )
                {
                    it->second[TypeID::L3] -= satPCenterL3;
                }

                // Look for C5
                if( it->second.find(TypeID::C5) != it->second.end() )
                {
                    it->second[TypeID::C5] -= satPCenterL5;
                }

                // Look for L5
                if( it->second.find(TypeID::L5) != it->second.end() )
                {
                    it->second[TypeID::L5] -= satPCenterL5;
                }

                // Look for C6
                if( it->second.find(TypeID::C6) != it->second.end() )
                {
                    it->second[TypeID::C6] -= satPCenterL6;
                }

                // Look for L6
                if( it->second.find(TypeID::L6) != it->second.end() )
                {
                    it->second[TypeID::L6] -= satPCenterL6;
                }

                // Look for C7
                if( it->second.find(TypeID::C7) != it->second.end() )
                {
                    it->second[TypeID::C7] -= satPCenterL7;
                }

                // Look for L7
                if( it->second.find(TypeID::L7) != it->second.end() )
                {
                    it->second[TypeID::L7] -= satPCenterL7;
                }

                // Look for C8
                if( it->second.find(TypeID::C8) != it->second.end() )
                {
                    it->second[TypeID::C8] -= satPCenterL8;
                }

                // Look for L8
                if( it->second.find(TypeID::L8) != it->second.end() )
                {
                    it->second[TypeID::L8] -= satPCenterL8;
                }

                // Look for C9
                if( it->second.find(TypeID::C9) != it->second.end() )
                {
                    it->second[TypeID::C9] -= satPCenterL9;
                }

                // Look for L9
                if( it->second.find(TypeID::L9) != it->second.end() )
                {
                    it->second[TypeID::L9] -= satPCenterL9;
                }

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

    }  // End of method 'ComputeSatPCenter::Process()'


     /** Return a gnssDataMap object, adding the new data generated
      *  when calling this object.
      *
      * @param gData    Data object holding the data.
      */
    gnssDataMap& ComputeSatPCenter::Process(gnssDataMap& gData)
        throw(ProcessingException)
    {
        SourceIDSet sourceRejectedSet;

        if(pEphStore == NULL)
        {
            ProcessingException e("Pointer to EphemerisStore NULL!");
            GPSTK_THROW(e);
        }

        if(pAntexReader == NULL)
        {
            ProcessingException e("Pointer to AntexReader NULL!");
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

                Process( gdmIt->first, sdmIt->second );
//                break;
            }
//            break;
        }

        gData.removeSourceID( sourceRejectedSet );

        return gData;

    }  // End of method 'ComputeSatPCenter::Process()'


}  // End of namespace gpstk
