//============================================================================
//
//  This file is part of GPSTk, the GPS Toolkit.
//
//  The GPSTk is free software; you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published
//  by the Free Software Foundation; either version 3.0 of the License, or
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
//  Copyright 2004, The University of Texas at Austin
//
//  Kaifa Kuang - Wuhan University . 2015
//
//============================================================================

//============================================================================
//
//This software developed by Applied Research Laboratories at the University of
//Texas at Austin, under contract to an agency or agencies within the U.S.
//Department of Defense. The U.S. Government retains all rights to use,
//duplicate, distribute, disclose, or release this software.
//
//Pursuant to DoD Directive 523024
//
// DISTRIBUTION STATEMENT A: This software has been approved for public
//                           release, distribution is unlimited.
//
//=============================================================================

/**
 * @file EGM08Model.cpp
 */

#include "EGM08Model.hpp"
#include "constants.hpp"
#include "StringUtils.hpp"
#include "Legendre.hpp"
#include "Epoch.hpp"


using namespace std;
using namespace gpstk::StringUtils;


namespace gpstk
{
    /// Load file
    void EGM08Model::loadFile(const string& file)
        throw(FileMissingException)
    {
        ifstream inpf(file.c_str());

        if(!inpf)
        {
            FileMissingException fme("Could not open EGM file " + file);
            GPSTK_THROW(fme);
        }

        // First, file header
        string temp;
        while( getline(inpf,temp) )
        {
            if(temp.substr(0,11) == "end_of_head") break;
        }

        bool ok(true);

        string line;

        // Then, file data
        while( !inpf.eof() && inpf.good() )
        {
            getline(inpf,line);
            stripTrailing(line,'\r');

            if( inpf.eof() ) break;

            if( inpf.bad() ) { ok = false; break; }

            // degree, order
            int L, M;
            L        =  asInt( line.substr( 5, 4) );
            M        =  asInt( line.substr(10, 4) );

            // Cnm, Snm, sigmaCnm, sigmaSnm
            double C, S, sigmaC, sigmaS;
            C        =  for2doub( line.substr(17,22) );
            S        =  for2doub( line.substr(42,22) );
            sigmaC   =  for2doub( line.substr(68,16) );
            sigmaS   =  for2doub( line.substr(88,16) );

            int id = indexTranslator(L,M)-1;

//            cout << line << endl;

            if(L<=desiredDegree && M<=desiredOrder)
            {
                egmData.normalizedCS(id, 0) = C;
                egmData.normalizedCS(id, 1) = S;
                egmData.normalizedCS(id, 2) = sigmaC;
                egmData.normalizedCS(id, 3) = sigmaS;
            }
            else
            {
                break;
            }

        }  // End of 'while(...)'

        inpf.close();

        if( !ok )
        {
            FileMissingException fme("EGM file " + file + " is corrupted or in wrong format");
            GPSTK_THROW(fme);
        }

    }  // End of method 'EGM08Model::loadFile()'


    /** Compute acceleration (and related partial derivatives) of EGM.
     * @param tt        TT
     * @param orbits    orbits
     */
    void EGM08Model::Compute( const CommonTime&   tt,
                              const satVectorMap& orbits )
    {

        // make a copy of Spherical Harmonic Coefficients
        Matrix<double> CS( egmData.normalizedCS );

        // compute time in years since J2000
        CommonTime utc( pRefSys->TT2UTC(tt) );
        double MJD_UTC = MJD(utc).mjd;
        double ly1 = (MJD_UTC - 51544.0) / 365.25;
        double ly2 = ly1 * ly1;
        double ly3 = ly2 * ly1;

        // Low-degree coefficients of the conventional geopotential model
        // see IERS Conventions 2010, Table 6.2
        // Note: the term C20 has been replaced with the corresponding
        // tide-free one
        const double value[3] =
        {
            // the value of C20, C30 and C40 at 2000.0
            -0.48416531e-3, 0.9571612e-6, 0.5399659e-6
        };

        const double rate[3]  =
        {
            // the rate of C20, C30 and C40
            11.6e-12,       4.9e-12,      4.7e-12
        };


        // indexes for degree = 2
        int id20 = indexTranslator(2,0) - 1;
        int id21 = indexTranslator(2,1) - 1;
        int id22 = indexTranslator(2,2) - 1;

        // indexes for degree = 3
        int id30 = indexTranslator(3,0) - 1;
        int id31 = indexTranslator(3,1) - 1;
        int id32 = indexTranslator(3,2) - 1;
        int id33 = indexTranslator(3,3) - 1;

        // indexes for degree = 4
        int id40 = indexTranslator(4,0) - 1;
        int id41 = indexTranslator(4,1) - 1;
        int id42 = indexTranslator(4,2) - 1;
        int id43 = indexTranslator(4,3) - 1;
        int id44 = indexTranslator(4,4) - 1;


        // The instantaneous value of coefficients Cn0 to be used when computing
        // orbits
        // see IERS Conventions 2010, Equations 6.4
        CS(id20, 0) = value[0] + ly1*rate[0];  // C20
        CS(id30, 0) = value[1] + ly1*rate[1];  // C30
        CS(id40, 0) = value[2] + ly1*rate[2];  // C40


        // Coefficients of the IERS (2010) mean pole model
        // see IERS Conventions 2010, Table 7.7

        // until 2010.0, unit: mas/year
        const double xp1[4] = {  55.974, 1.8243,  0.18413,  0.007024 };
        const double yp1[4] = { 346.346, 1.7896, -0.10729, -0.000908 };
        // after 2010.0, unit: mas/year
        const double xp2[2] = {  23.513,  7.6141 };
        const double yp2[2] = { 358.891, -0.6287 };

        // get the mean pole at epoch 2000.0 from IERS Conventions 2010
        // see IERS Conventions 2010, Equation 7.25
        double xpm(0.0), ypm(0.0);

        if(MJD_UTC < 55197.0)    // until 2010.0
        {
            xpm = ( xp1[0] + xp1[1]*ly1 + xp1[2]*ly2 + xp1[3]*ly3 )*1e-3;
            ypm = ( yp1[0] + yp1[1]*ly1 + yp1[2]*ly2 + yp1[3]*ly3 )*1e-3;
        }
        else                     // after 2010.0
        {
            xpm = ( xp2[0] + xp2[1]*ly1 )*1e-3;
            ypm = ( yp2[0] + yp2[1]*ly1 )*1e-3;
        }

        // convert pole position from arcseconds to radians
        xpm = xpm * AS_TO_RAD;
        ypm = ypm * AS_TO_RAD;

        // Rotate from the Earth-fixed frame, where the coefficients are
        // pertinent, to an inertial frame, where the satellite motion is
        // computed
        // see IERS Conventions 2010, Equation 6.5
        //
        // C21 = +sqrt(3)*xpm*C20 - xpm*C22 + ypm*S22
        // S21 = -sqrt(3)*ypm*C20 - ypm*C22 - xpm*S22
        //
        double C20 = CS(id20,0);
        double C22 = CS(id22,0); double S22 = CS(id22,1);

        double C21 = +std::sqrt(3.0)*xpm*C20 - xpm*C22 + ypm*S22;
        double S21 = -std::sqrt(3.0)*ypm*C20 - ypm*C22 - xpm*S22;

        CS(id21,0) = C21; CS(id21,1) = S21;


        //// Tide corrections ////

        Matrix<double> dCS;

        // solid Earth tides
        if(pSolidTide != NULL)
        {
            // corrections of CS
            dCS = pSolidTide->getSolidTide(tt);

            for(int i=0; i<dCS.rows(); ++i)
            {
                CS(i,0) += dCS(i,0);
                CS(i,1) += dCS(i,1);
            }
        }

        // ocean tides
        if(pOceanTide != NULL)
        {
            dCS = pOceanTide->getOceanTide(tt);

            for(int i=0; i<dCS.rows(); ++i)
            {
                CS(i,0) += dCS(i,0);
                CS(i,1) += dCS(i,1);
            }
        }

        // solid Earth pole tide and ocean pole tide
        if(pPoleTide != NULL)
        {
            dCS = pPoleTide->getPoleTide(tt);

            for(int i=0; i<dCS.rows(); ++i)
            {
                CS(i,0) += dCS(i,0);
                CS(i,1) += dCS(i,1);
            }
        }


        // transformation matrixes between ECI and ECF
        Matrix<double> C2T( pRefSys->C2TMatrix(utc) );
        Matrix<double> T2C( transpose(C2T) );

        SatID sat;
        Vector<double> orbit;

        Vector<double> r_sat_eci(3,0.0);
        Vector<double> r_sat_ecf(3,0.0);
        double rx(0.0), ry(0.0), rz(0.0);
        double rho(0.0), lat(0.0), lon(0.0);
        double slat(0.0), clat(0.0);

        Vector<double> slon(desiredDegree+1,0.0);
        Vector<double> clon(desiredDegree+1,0.0);
        slon(0) = 0.0; clon(0) = 1.0;

        Vector<double> leg0, leg1, leg2;

        // partials of (rho, lat, lon) in ECF to (x, y, z) in ECF
        Matrix<double> b(3,3,0.0);

        Matrix<double> db_drho(3,3,0.0);     // db / drho
        Matrix<double> db_dlat(3,3,0.0);     // db / dlat
        Matrix<double> db_dlon(3,3,0.0);     // db / dlon

        ///////////// partials of v to (rho, lat, lon) in ECF /////////////////
        //                                                                   //
        //                           | dv / drho |                           //
        //                  f_rll =  | dv / dlat |                           //
        //                           | dv / dlon |                           //
        //                                                                   //
        ///////////////////////////////////////////////////////////////////////
        Vector<double> f_rll(3,0.0);


        ///////////// partials of f_rll to (rho, lat, lon) in ECF /////////////
        //                                                                   //
        //         | df_rll(0) / drho, df_rll(0) / dlat, df_rll(0) / dlon |  //
        // g_rll = | df_rll(1) / drho, df_rll(1) / dlat, df_rll(1) / dlon |  //
        //         | df_rll(2) / drho, df_rll(2) / dlat, df_rll(2) / dlon |  //
        //                                                                   //
        ///////////////////////////////////////////////////////////////////////
        Matrix<double> g_rll(3,3,0.0);


        double gm_r1(0.0), gm_r2(0.0), gm_r3(0.0);

        // (ae/rho)^n
        double fct(0.0);

        // index for (n,m)
        int idnm(0);

        // Pnm and its derivatives wrt latitude
        double P0(0.0), P1(0.0), P2(0.0);

        // Cnm, Snm
        double Cnm(0.0), Snm(0.0);

        // sin(m*lon) and cos(m*lon)
        double smlon(0.0), cmlon(0.0);

        // gravitation acceleration in ECI
        Vector<double> a(3,0.0);

        // partials of gravitation acceleration in ECI to (x, y, z) in ECI
        Matrix<double> df_ecf_drll(3,3,0.0);

        Matrix<double> da_dr(3,3,0.0);

        for( satVectorMap::const_iterator it = orbits.begin();
             it != orbits.end();
             ++it )
        {
            sat = it->first;
            orbit = it->second;

            // satellite position in ECI
            r_sat_eci(0) = orbit(0);
            r_sat_eci(1) = orbit(1);
            r_sat_eci(2) = orbit(2);

            // satellite position in ECF
            r_sat_ecf = C2T * r_sat_eci;
            rx = r_sat_ecf(0);
            ry = r_sat_ecf(1);
            rz = r_sat_ecf(2);

            // geocentric distance, latitude and longitude of satellite
            // the latitude is in [-PI/2,+PI/2], the longitude is in [-PI,+PI]
            rho = norm(r_sat_ecf);
            lat = std::atan( rz / std::sqrt(rx*rx + ry*ry) );
            lon = std::atan2( ry, rx );

            // sine and cosine of geocentric latitude
            slat = std::sin(lat);
            clat = std::cos(lat);

            // sine and cosine of geocentric longitude
            //    slon(0) = sin(0*lon);
            //    slon(1) = sin(1*lon);
            //        ......
            //    slon(i) = slon(i-1)*cos(lon) + clon(i-1)*sin(lon)
            //    clon(i) = clon(i-1)*cos(lon) - slon(i-1)*sin(lon)
            slon(1) = std::sin(lon);  clon(1) = std::cos(lon);

            for(int i=2; i<=desiredDegree; ++i)
            {
                slon(i) = slon(i-1)*clon(1) + clon(i-1)*slon(1);
                clon(i) = clon(i-1)*clon(1) - slon(i-1)*slon(1);
            }

            // fully normalized associated legendre functions and its gradients
            legendre(desiredDegree, lat, leg0, leg1, leg2, 2);

            // partials of (rho, lat, lon) in ITRS to (x, y, z) in ECF
            b(0,0) =  clat * clon(1);
            b(0,1) =  clat * slon(1);
            b(0,2) =  slat;
            b(1,0) = -slat * clon(1) / rho;
            b(1,1) = -slat * slon(1) / rho;
            b(1,2) =  clat / rho;
            b(2,0) = -slon(1) / (rho * clat);
            b(2,1) =  clon(1) / (rho * clat);


            // partials of b to (rho, lat, lon) in ECF
            db_drho(1,0) =  slat * clon(1) / (rho*rho);
            db_drho(1,1) =  slat * slon(1) / (rho*rho);
            db_drho(1,2) = -clat / (rho*rho);
            db_drho(2,0) =  slon(1) / (rho*rho * clat);
            db_drho(2,1) = -clon(1) / (rho*rho * clat);

            db_dlat(0,0) = -slat * clon(1);
            db_dlat(0,1) = -slat * slon(1);
            db_dlat(0,2) =  clat;
            db_dlat(1,0) = -clat * clon(1) / rho;
            db_dlat(1,1) = -clat * slon(1) / rho;
            db_dlat(1,2) = -slat / rho;
            db_dlat(2,0) = -slat * slon(1) / (rho * clat*clat);
            db_dlat(2,1) =  slat * clon(1) / (rho * clat*clat);

            db_dlon(0,0) = -clat * slon(1);
            db_dlon(0,1) =  clat * clon(1);
            db_dlon(1,0) =  slat * slon(1) / rho;
            db_dlon(1,1) = -slat * clon(1) / rho;
            db_dlon(2,0) = -clon(1) / (rho * clat);
            db_dlon(2,1) = -slon(1) / (rho * clat);


            // (GM/rho)^1, (GM/rho)^2, (GM/rho)^3
            gm_r1 = egmData.GM / rho;
            gm_r2 = gm_r1 / rho;
            gm_r3 = gm_r2 / rho;

            f_rll.resize(3,0.0);
            g_rll.resize(3,3,0.0);

            // loop for degree
            for(int n=0; n<=desiredDegree; ++n)
            {
                fct = std::pow(egmData.ae/rho, n);

                // loop for order
                for(int m=0; m<=n && m<=desiredOrder; ++m)
                {
                    idnm = indexTranslator(n,m) - 1;

                    P0 = leg0(idnm);
                    P1 = leg1(idnm);
                    P2 = leg2(idnm);

                    Cnm = CS(idnm,0);
                    Snm = CS(idnm,1);

                    smlon = slon(m);
                    cmlon = clon(m);

                    // f_rll
                    f_rll(0) += -gm_r2*(n+1)*fct*P0*(Cnm*cmlon+Snm*smlon);
                    f_rll(1) +=  gm_r1*fct*P1*(Cnm*cmlon+Snm*smlon);
                    f_rll(2) +=  gm_r1*m*fct*P0*(-Cnm*smlon+Snm*cmlon);

                    // g_rll
                    g_rll(0,0) +=  gm_r3*(n+1)*(n+2)*fct*P0*(Cnm*cmlon+Snm*smlon);
                    g_rll(0,1) += -gm_r2*(n+1)*fct*P1*(Cnm*cmlon+Snm*smlon);
                    g_rll(0,2) += -gm_r2*(n+1)*m*fct*P0*(-Cnm*smlon+Snm*cmlon);
                    g_rll(1,1) +=  gm_r1*fct*P2*(Cnm*cmlon+Snm*smlon);
                    g_rll(1,2) +=  gm_r1*m*fct*P1*(-Cnm*smlon+Snm*cmlon);
                    g_rll(2,2) += -gm_r1*m*m*fct*P0*(Cnm*cmlon+Snm*smlon);

                }   // End of 'for(int m=0; ...)'

            }   // End of 'for(int n=0; ...)'


            // symmetry of gradiometry
            g_rll(1,0) = g_rll(0,1);
            g_rll(2,0) = g_rll(0,2);
            g_rll(2,1) = g_rll(1,2);

            a = T2C * transpose(b) * f_rll;
            satAcc[sat] = a;

            for(int i=0; i<3; ++i)
            {
                df_ecf_drll(i,0) = g_rll(0,0)*b(0,i) + f_rll(0)*db_drho(0,i)
                                 + g_rll(1,0)*b(1,i) + f_rll(1)*db_drho(1,i)
                                 + g_rll(2,0)*b(2,i) + f_rll(2)*db_drho(2,i);
                df_ecf_drll(i,1) = g_rll(0,1)*b(0,i) + f_rll(0)*db_dlat(0,i)
                                 + g_rll(1,1)*b(1,i) + f_rll(1)*db_dlat(1,i)
                                 + g_rll(2,1)*b(2,i) + f_rll(2)*db_dlat(2,i);
                df_ecf_drll(i,2) = g_rll(0,2)*b(0,i) + f_rll(0)*db_dlon(0,i)
                                 + g_rll(1,2)*b(1,i) + f_rll(1)*db_dlon(1,i)
                                 + g_rll(2,2)*b(2,i) + f_rll(2)*db_dlon(2,i);
            }

            da_dr = T2C * df_ecf_drll * b * C2T;
            satPartialR[sat] = da_dr;

        } // End of 'for(satVectorMap::const_iterator...)'

    } // End of method 'EGM08Model::Compute(...)'

}   // End of namespace 'gpstk'
