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

/**
 * @file SRPModel.cpp
 * Class to do SRP calculation.
 */

#include "SRPModel.hpp"
#include "constants.hpp"

namespace gpstk
{

    /* Determines if the satellite is in sunlight or shadow.
     * Taken from Montenbruck and Gill p. 80-83
     * @param r_sat     position of satellite (wrt Earth) [m]
     * @param r_sun     position of Sun (wrt Earth) [m]
     * @param r_moon    position of Moon (wrt Earth) [m]
     * @return          0.0 if in shadow
     *                  1.0 if in sunlight
     *                  0 to 1.0 if in partial shadow
     */
    double SRPModel::getShadowFunction(Vector<double> r_sat,
                                       Vector<double> r_sun,
                                       Vector<double> r_moon,
                                       SRPModel::ShadowModel sm)
    {
        // shadow function
        double v( 0.0 );

        // Sun direction unit vector
        Vector<double> e_sun( normalize(r_sun) );

        double r_dot_sun( dot(r_sat,e_sun) );

        if(r_dot_sun>0)
        {
            // Sunny side of central body is always fully lit and return
            v= 1.0;
            return v;
        }

        if(sm == SM_CYLINDRICAL)
        {
            // Taken fram Jisheng Li P111, and checked with GMAT and Bernese5 SHADOW.f
            v = ((r_dot_sun>0 || norm(r_sat-r_dot_sun*e_sun)>RE_EARTH) ? 1.0 : 0.0);
            return v;
        }
        else if(sm == SM_CONICAL)
        {
            double r_mag( norm(r_sat) );

            // vector from satellite to Sun
            Vector<double> d( r_sun - r_sat );

            double dmag = norm(d);

            double a = std::asin(RE_SUN/dmag);
            double b = std::asin(RE_EARTH/r_mag);
            double c = std::acos(-1.0*dot(r_sat,d)/(r_mag*dmag));

            if((a+b) <= c)         // in Sun light
            {
                v = 1.0;
            }
            else if(c < (b-a))     // in Umbra
            {
                v =  0.0;
            }
            else                   // in Penumbra
            {
                double x( (c*c+a*a-b*b)/(2*c) );
                double y( std::sqrt(a*a-x*x) );
                double A( a*a*std::acos(x/a)+b*b*std::acos((c-x)/b)-c*y );
                v = 1.0 - A/(PI*a*a);
            }

            return v;
        }
        else
        {
            // unexpected value
            Exception e("Unexpect ShadowModel in getShadowFunction()");
            GPSTK_THROW(e);
        }

        return v;

    }  // End of method 'SRPModel::getShadowFunction()'


}  // End of namespace 'gpstk'
