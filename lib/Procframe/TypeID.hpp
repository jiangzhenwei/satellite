#pragma ident "$Id$"

/**
 * @file TypeID.hpp
 * gpstk::TypeID - This class was written taking as inspiration ObsID. The
 * objective of this class is to create an index able to represent any type
 * of observation, correction, model parameter or other data value of interest
 * for GNSS data processing. This class is extensible in run-time, so the
 * programmer may add indexes on-demand.
 */

#ifndef GPSTK_TYPEID_HPP
#define GPSTK_TYPEID_HPP

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
//  Dagoberto Salazar - gAGE ( http://www.gage.es ). 2007, 2008, 2009
//
//============================================================================


#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <map>
#include "RinexObsHeader.hpp"
#include "Rinex3ObsHeader.hpp"
#include "RinexObsID.hpp"

namespace gpstk
{

    /** This class creates an index able to represent any type of observation,
     *  correction, model parameter or other data value of interest for GNSS
     *  data processing.
     *
     * This class is extensible in run-time, so the programmer may add
     * indexes on-demand. For instance, in order to create a new TypeID
     * object referring INS-related data, and with "Inertial" as description
     * string, you may write the following:
     *
     * @code
     *    TypeID INS = TypeID::newValueType("Inertial");
     * @endcode
     *
     * Or using the constructor:
     *
     * @code
     *    TypeID INS(TypeID::newValueType("Inertial"));
     * @endcode
     *
     * From now on, you'll be able to use INS as TypeID when you need to
     * refer to inertial system data.
     *
     */
    class TypeID
    {
    public:

        /// The type of the data value.
        enum ValueType
        {
            Unknown,

            // C1*
            C1C,
            C1S,
            C1L,
            C1X,
            C1P,
            C1W,
            C1Y,
            C1M,
            C1A,
            C1B,
            C1Z,

            // L1*
            L1C,
            L1S,
            L1L,
            L1X,
            L1P,
            L1W,
            L1Y,
            L1M,
            L1N,
            L1A,
            L1B,
            L1Z,

            // D1*
            D1C,
            D1S,
            D1L,
            D1X,
            D1P,
            D1W,
            D1Y,
            D1M,
            D1N,
            D1A,
            D1B,
            D1Z,

            // S1*
            S1C,
            S1S,
            S1L,
            S1X,
            S1P,
            S1W,
            S1Y,
            S1M,
            S1N,
            S1A,
            S1B,
            S1Z,

            // C2*
            C2C,
            C2D,
            C2S,
            C2L,
            C2X,
            C2P,
            C2W,
            C2Y,
            C2M,
            C2I,
            C2Q,

            // L2*
            L2C,
            L2D,
            L2S,
            L2L,
            L2X,
            L2P,
            L2W,
            L2Y,
            L2M,
            L2N,
            L2I,
            L2Q,

            // D2*
            D2C,
            D2D,
            D2S,
            D2L,
            D2X,
            D2P,
            D2W,
            D2Y,
            D2M,
            D2N,
            D2I,
            D2Q,

            // S2*
            S2C,
            S2D,
            S2S,
            S2L,
            S2X,
            S2P,
            S2W,
            S2Y,
            S2M,
            S2N,
            S2I,
            S2Q,

            // C3*
            C3I,
            C3Q,
            C3X,

            // L3*
            L3I,
            L3Q,
            L3X,

            // D3*
            D3I,
            D3Q,
            D3X,

            // S3*
            S3I,
            S3Q,
            S3X,

            // C5*
            C5I,
            C5Q,
            C5X,
            C5A,
            C5B,
            C5C,

            // L5*
            L5I,
            L5Q,
            L5X,
            L5A,
            L5B,
            L5C,

            // D5*
            D5I,
            D5Q,
            D5X,
            D5A,
            D5B,
            D5C,

            // S5*
            S5I,
            S5Q,
            S5X,
            S5A,
            S5B,
            S5C,

            // C6*
            C6A,
            C6B,
            C6C,
            C6X,
            C6Z,
            C6S,
            C6L,
            C6I,
            C6Q,

            // L6*
            L6A,
            L6B,
            L6C,
            L6X,
            L6Z,
            L6S,
            L6L,
            L6I,
            L6Q,

            // D6*
            D6A,
            D6B,
            D6C,
            D6X,
            D6Z,
            D6S,
            D6L,
            D6I,
            D6Q,

            // S6*
            S6A,
            S6B,
            S6C,
            S6X,
            S6Z,
            S6S,
            S6L,
            S6I,
            S6Q,

            // C7*
            C7I,
            C7Q,
            C7X,

            // L7*
            L7I,
            L7Q,
            L7X,

            // D7*
            D7I,
            D7Q,
            D7X,

            // S7*
            S7I,
            S7Q,
            S7X,

            // C8*
            C8I,
            C8Q,
            C8X,

            // L8*
            L8I,
            L8Q,
            L8X,

            // D8*
            D8I,
            D8Q,
            D8X,

            // S8*
            S8I,
            S8Q,
            S8X,

            // C9*
            C9A,
            C9B,
            C9C,
            C9X,

            // L9*
            L9A,
            L9B,
            L9C,
            L9X,

            // D9*
            D9A,
            D9B,
            D9C,
            D9X,

            // S9*
            S9A,
            S9B,
            S9C,
            S9X,


            // LLI
            LLI1,
            LLI2,
            LLI3,
            LLI5,
            LLI6,
            LLI7,
            LLI8,
            LLI9,

            // SSI
            SSI1,
            SSI2,
            SSI3,
            SSI5,
            SSI6,
            SSI7,
            SSI8,
            SSI9,


            // C*
            C1,
            C2,
            C3,
            C5,
            C6,
            C7,
            C8,
            C9,
            P1,
            P2,
            P3,
            P5,
            P6,
            P7,
            P8,
            P9,

            // L*
            L1,
            L2,
            L3,
            L4,
            L5,
            L6,
            L7,
            L8,
            L9,

            // D*
            D1,
            D2,
            D3,
            D5,
            D6,
            D7,
            D8,
            D9,

            // S*
            S1,
            S2,
            S3,
            S5,
            S6,
            S7,
            S8,
            S9,


            // RINEX 2.20
            LA,     // Phase measurements on L1 derived from C/A code tracking
            DA,     // Doppler frequency on LA
            T1,     // Transit Integrated Doppler on 150 (T1) and 400 MHz (T2)
            T2,
            SA,     // Raw signal or SNR for LA
            CH,     // Receiver channel number


            // Cycle Slip
            CSL1,
            CSL2,
            CSL3,
            CSL5,
            CSL6,
            CSL7,
            CSL8,
            CSL9,

            // Linear Combination
            PC,             ///< Ionospheric-free code combination
            LC,             ///< Ionospheric-free phase combination
            PI,             ///< Ionospheric code combination
            LI,             ///< Ionospheric phase combination
            PW,             ///< Wide-lane code combination
            LW,             ///< Wide-lane phase combination
            PN,             ///< Narrow-lane code combination
            LN,             ///< Narrow-lane phase combination
            MW,             ///< Melbourne-Wubbena combination
            WL,             ///< Wide-lane combination(+1*L1-1*L2)
            WL1,            ///< Wide-lane combination(-1*L1+2*L2)
            WL2,            ///< Wide-lane combination(-2*L1+3*L2)
            WL3,            ///< Wide-lane combination(-3*L1+4*L2)
            WL4,            ///< Wide-lane combination(+4*L1-5*L2)
            EWL,            ///< Wide-lane combination(-7*L1+9*L2)

            Q1,             ///< precise code with minus ionospheric delays in L1
            Q2,             ///< precise code with minus ionospheric delays in L2
            Q3,             ///< precise code with minus ionospheric delays in L3
            Q5,             ///< precise code with minus ionospheric delays in L5
            Q6,             ///< precise code with minus ionospheric delays in L6
            Q7,             ///< precise code with minus ionospheric delays in L7
            Q8,             ///< precise code with minus ionospheric delays in L8
            Q9,             ///< precise code with minus ionospheric delays in L9

            PC12,           ///< Ionospheric-free code combination, L1+L2
            LC12,           ///< Ionospheric-free phase combination, L1+L2
            PC15,           ///< Ionospheric-free code combination, L1+L5
            LC15,           ///< Ionospheric-free phase combination, L1+L5
            PC27,           ///< Ionospheric-free code combination, L2+L7
            LC27,           ///< Ionospheric-free phase combination, L2+L7

            PI12,           ///< Ionospheric code combination, L1+L2
            LI12,           ///< Ionospheric phase combination, L1+L2
            PI15,           ///< Ionospheric code combination, L1+L5
            LI15,           ///< Ionospheric phase combination, L1+L5
            PI27,           ///< Ionospheric code combination, L2+L7
            LI27,           ///< Ionospheric phase combination, L2+L7

            PW12,           ///< Wide-lane code combination, L1+L2
            LW12,           ///< Wide-lane phase combination, L1+L2
            PW15,           ///< Wide-lane code combination, L1+L5
            LW15,           ///< Wide-lane phase combination, L1+L5
            PW27,           ///< Wide-lane code combination, L2+L7
            LW27,           ///< Wide-lane phase combination, L2+L7

            PN12,           ///< Narrow-lane code combination, L1+L2
            LN12,           ///< Narrow-lane phase combination, L1+L2
            PN15,           ///< Narrow-lane code combination, L1+L5
            LN15,           ///< Narrow-lane phase combination, L1+L5
            PN27,           ///< Narrow-lane code combination, L2+L7
            LN27,           ///< Narrow-lane phase combination, L2+L7

            MW12,           ///< Melbourne-Wubbena combination, L1+L2
            MW15,           ///< Melbourne-Wubbena combination, L1+L5
            MW27,           ///< Melbourne-Wubbena combination, L2+L7

            Q1OfL1L2,
            Q2OfL1L2,

            Q1OfL1L5,
            Q5OfL1L5,

            Q2OfL2L7,
            Q7OfL2L7,


            // Prefit Residual
            prefitC1,
            prefitC2,
            prefitC3,
            prefitC5,
            prefitC6,
            prefitC7,
            prefitC8,
            prefitC9,
            prefitC,
            prefitCWithoutClock,
            prefitCWithSatClock,
            prefitCWithStaClock,
            prefitCForPCE,
            prefitCForPOD,
            prefitCForPPP,


            prefitL1,
            prefitL2,
            prefitL3,
            prefitL5,
            prefitL6,
            prefitL7,
            prefitL8,
            prefitL9,
            prefitL,
            prefitLWithoutClock,
            prefitLWithSatClock,
            prefitLWithStaClock,
            prefitLForPCE,
            prefitLForPOD,
            prefitLForPPP,

            prefitPC,
            prefitLC,
            prefitMW,
            prefitWL,
            prefitWL1,
            prefitWL2,
            prefitWL3,
            prefitWL4,
            prefitEWL,


            prefitC12,
            prefitC12WithoutClock,
            prefitC12WithSatClock,
            prefitC12WithStaClock,
            prefitC12ForPCE,
            prefitC12ForPOD,
            prefitC12ForPPP,

            prefitC15,
            prefitC15WithoutClock,
            prefitC15WithSatClock,
            prefitC15WithStaClock,
            prefitC15ForPCE,
            prefitC15ForPOD,
            prefitC15ForPPP,

            prefitC27,
            prefitC27WithoutClock,
            prefitC27WithSatClock,
            prefitC27WithStaClock,
            prefitC27ForPCE,
            prefitC27ForPOD,
            prefitC27ForPPP,


            prefitL12,
            prefitL12WithoutClock,
            prefitL12WithSatClock,
            prefitL12WithStaClock,
            prefitL12ForPCE,
            prefitL12ForPOD,
            prefitL12ForPPP,

            prefitL15,
            prefitL15WithoutClock,
            prefitL15WithSatClock,
            prefitL15WithStaClock,
            prefitL15ForPCE,
            prefitL15ForPOD,
            prefitL15ForPPP,

            prefitL27,
            prefitL27WithoutClock,
            prefitL27WithSatClock,
            prefitL27WithStaClock,
            prefitL27ForPCE,
            prefitL27ForPOD,
            prefitL27ForPPP,

            prefitMW12,
            prefitMW15,
            prefitMW27,


            // Postfit Residual
            postfitC1,
            postfitC2,
            postfitC3,
            postfitC5,
            postfitC6,
            postfitC7,
            postfitC8,
            postfitC9,
            postfitC,

            postfitL1,
            postfitL2,
            postfitL3,
            postfitL5,
            postfitL6,
            postfitL7,
            postfitL8,
            postfitL9,
            postfitL,

            postfitPC,
            postfitLC,
            postfitMW,
            postfitWL,
            postfitWL1,
            postfitWL2,
            postfitWL3,
            postfitWL4,
            postfitEWL,

            postfitC12,
            postfitC15,
            postfitL12,
            postfitL15,

            // Phase Ambiguity
            BL1,
            BL2,
            BL3,
            BL5,
            BL6,
            BL7,
            BL8,
            BL9,
            BLC,
            BWL,
            BWL1,
            BWL2,
            BWL3,
            BWL4,

            BLC12,
            BLC15,

            BWL12,
            BWL15,


            // Uncalibrated Phase Delay
            updStaL1,
            updStaL2,
            updStaL3,
            updStaL5,
            updStaL6,
            updStaL7,
            updStaL8,
            updStaL9,
            updStaLC12,
            updStaWL12,

            updSatL1,
            updSatL2,
            updSatL3,
            updSatL5,
            updSatL6,
            updSatL7,
            updSatL8,
            updSatL9,
            updSatLC12,
            updSatWL12,

            // Instrumental Delay
            instC1,
            instC2,
            instC3,
            instC5,
            instC6,
            instC7,
            instC8,
            instC9,
            instPC,
            instPN,

            instL1,
            instL2,
            instL3,
            instL5,
            instL6,
            instL7,
            instL8,
            instL9,
            instLC,
            instLW,

            instMW,


            instPC12,
            instPC15,
            instPC27,
            instLC12,
            instLC15,
            instLC27,

            instPN12,
            instPN15,
            instPN27,
            instLW12,
            instLW15,
            instLW27,

            instMW12,
            instMW15,
            instMW27,


            meanMW12,
            varMW12,


            // ISB
            ISB_GLO,
            ISB_Gal,
            ISB_BDS,


            // Geometric Distance
            rho,
            rhoDot,
            rhoDot2,


            // Relativity
            relativity,


            // Gravitational Delay
            gravDelay,


            // Tropospheric Delay
            tropo,          ///< Vertical tropospheric delay, total
            tropoWeight,    ///< Vertical tropospheric delay weight, total
            dryTropo,       ///< Vertical tropospheric delay, dry component
            dryMap,         ///< Tropospheric mapping function, dry component
            wetTropo,       ///< Vertical tropospheric delay, wet component
            wetMap,         ///< Tropospheric mapping function, wet component
            tropoSlant,     ///< Slant tropospheric delay, total


            // Ionospheric Delay
            iono,           ///< Vertical ionospheric delay
            ionoTEC,        ///< Total Electron Content (in TECU), 1TECU = 1e+16 electrons per m**2
            ionoMap,        ///< Ionospheric mapping function
            ionoMap2,       ///< Ionospheric mapping function for second order ionospheric delay
            ionoL1,         ///< Slant ionospheric delay, frequency L1
            ionoL2,         ///< Slant ionospheric delay, frequency L2
            ionoL3,         ///< Slant ionospheric delay, frequency L3
            ionoL5,         ///< Slant ionospheric delay, frequency L5
            ionoL6,         ///< Slant ionospheric delay, frequency L6
            ionoL7,         ///< Slant ionospheric delay, frequency L7
            ionoL8,         ///< Slant ionospheric delay, frequency L8
            ionoL9,         ///< Slant ionospheric delay, frequency L9
            ionoPN,         ///< Slant ionospheric delay, frequency PN
            ionoLW,         ///< Slant ionospheric delay, frequency LW
            ionoL1Weight,   ///< Weight for slant ionospheric delay on frequency L1


            // Epoch difference
            diffPrefitC1,
            diffPrefitC2,
            diffPrefitC,

            diffPrefitL1,
            diffPrefitL2,
            diffPrefitL,

            diffWetTropo,


            // Wind-Up Effect (in radians)
            windUp,


            // Eclipse Indicator
            eclipse,

            yawAngle,
            phiAngle,


            // Satellite Antenna Phase Center Correction
            satPCenter,     ///< Satellite antenna phase center correction, ray direction
            satPCenterX,    ///< Satellite antenna phase center correction, X component
            satPCenterY,    ///< Satellite antenna phase center correction, Y component
            satPCenterZ,    ///< Satellite antenna phase center correction, Z component


            // Satellite Antenna Phase Center Offset
            satPCOffX,      ///< Satellite antenna phase center offset, X component
            satPCOffY,      ///< Satellite antenna phase center offset, Y component
            satPCOffZ,      ///< Satellite antenna phase center offset, Z component


            // Station Antenna Phase Center Offset
            staPCOffU,      ///< Station antenna phase center offset, U component
            staPCOffN,      ///< Station antenna phase center offset, N component
            staPCOffE,      ///< Station antenna phase center offset, E component


            // Satellite Elevation and Azimuth
            elevation,
            azimuth,


            // Satellite Arc Number
            satArc,


            // Weight
            weight,         ///< Weight assigned to a given observation
            weightC,        ///< Weight assigned to a given code observation
            weightL,        ///< Weight assigned to a given phase observation


            // Station State
            staXECF,        ///< Station position in ECF, X component
            staYECF,        ///< Station position in ECF, Y component
            staZECF,        ///< Station position in ECF, Z component
            staVXECF,       ///< Station velocity in ECF, X component
            staVYECF,       ///< Station velocity in ECF, Y component
            staVZECF,       ///< Station velocity in ECF, Z component
            staAXECF,       ///< Station acceleration in ECF, X component
            staAYECF,       ///< Station acceleration in ECF, Y component
            staAZECF,       ///< Station acceleration in ECF, Z component

            staXECI,        ///< Station postion in ECI, X component
            staYECI,        ///< Station postion in ECI, Y component
            staZECI,        ///< Station postion in ECI, Z component
            staVXECI,       ///< Station velocity in ECI, X component
            staVYECI,       ///< Station velocity in ECI, Y component
            staVZECI,       ///< Station velocity in ECI, Z component
            staAXECI,       ///< Station acceleration in ECI, X component
            staAYECI,       ///< Station acceleration in ECI, Y component
            staAZECI,       ///< Station acceleration in ECI, Z component

            staLatECF,      ///< Station position, Latitude component
            staLonECF,      ///< Station position, Longitude component
            staHECF,        ///< Station position, Height component
            staVLatECF,     ///< Station velocity, Latitude component
            staVLonECF,     ///< Station velocity, Longitude component
            staVHECF,       ///< Station velocity, Height component
            staALatECF,     ///< Station acceleration, Latitude component
            staALonECF,     ///< Station acceleration, Longitude component
            staAHECF,       ///< Station acceleration, Height component

            // Station Clock
            cdtSta,         ///< Station clock offset
            cdtStaDot,      ///< Station clock drift
            cdtStaDot2,     ///< Station clock drift rate


            // Satellite State
            satXECF,        ///< Satellite position in ECF, X component
            satYECF,        ///< Satellite position in ECF, Y component
            satZECF,        ///< Satellite position in ECF, Z component
            satVXECF,       ///< Satellite velocity in ECF, X component
            satVYECF,       ///< Satellite velocity in ECF, Y component
            satVZECF,       ///< Satellite velocity in ECF, Z component
            satAXECF,       ///< Satellite acceleration in ECF, X component
            satAYECF,       ///< Satellite acceleration in ECF, Y component
            satAZECF,       ///< Satellite acceleration in ECF, Z component
            satXECI,        ///< Satellite position in ECI, X component
            satYECI,        ///< Satellite position in ECI, Y component
            satZECI,        ///< Satellite position in ECI, Z component
            satVXECI,       ///< Satellite velocity in ECI, X component
            satVYECI,       ///< Satellite velocity in ECI, Y component
            satVZECI,       ///< Satellite velocity in ECI, Z component
            satAXECI,       ///< Satellite acceleration in ECI, X component
            satAYECI,       ///< Satellite acceleration in ECI, Y component
            satAZECI,       ///< Satellite acceleration in ECI, Z component

            // Satellite Clock
            cdtSat,         ///< Satellite clock offset
            cdtSatDot,      ///< Satellite clock drift
            cdtSatDot2,     ///< Satellite clock drift rate
            cdtSatSum,      ///< Satellite clock offset sum


            // Correction or Coefficient of Station State
            dStaX,
            dStaY,
            dStaZ,
            dStaLat,
            dStaLon,
            dStaH,

            // Correction or Coefficient of Station Clock
            dcdtSta,
            dcdtStaDot,
            dcdtStaDot2,

            dcdtStaGPS,
            dcdtStaGal,

            ifcbSta,


            // Correction or Coefficient of Satellite State
            dSatX,
            dSatY,
            dSatZ,
            dSatVX,
            dSatVY,
            dSatVZ,

            dSatR,
            dSatT,
            dSatN,

            // Correction or Coefficient of Satellite Clock
            dcdtSat,
            dcdtSatDot,
            dcdtSatDot2,


            ifcbSat,


            // Solar Radiation Pressure (SRP) Coefficients
            SRPC1,
            SRPC2,
            SRPC3,
            SRPC4,
            SRPC5,
            SRPC6,
            SRPC7,
            SRPC8,
            SRPC9,

            // Correction or Coefficient of SRP Coefficients
            dSRPC1,
            dSRPC2,
            dSRPC3,
            dSRPC4,
            dSRPC5,
            dSRPC6,
            dSRPC7,
            dSRPC8,
            dSRPC9,


            // Earth Ratation Parameters (ERPs)
            xpole,          ///< X Pole
            ypole,          ///< Y Pole
            xpoleRate,      ///< X Pole Rate
            ypoleRate,      ///< Y Pole Rate
            UT1mUTC,        ///< UT1mUTC
            LOD,            ///< LOD

            // Correction or Coefficient of ERPs
            dxpole,
            dypole,
            dxpoleRate,
            dypoleRate,
            dUT1mUTC,
            dLOD,


            // Handy dummy types for non-standard processing
            dummy0,         ///< Generic, undefined type #0
            dummy1,         ///< Generic, undefined type #1
            dummy2,         ///< Generic, undefined type #2
            dummy3,         ///< Generic, undefined type #3
            dummy4,         ///< Generic, undefined type #4
            dummy5,         ///< Generic, undefined type #5
            dummy6,         ///< Generic, undefined type #6
            dummy7,         ///< Generic, undefined type #7
            dummy8,         ///< Generic, undefined type #8
            dummy9,         ///< Generic, undefined type #9

            Last,      ///< used to extend this...
            Placeholder = Last+1000
        };


        /// empty constructor, creates an invalid object
        TypeID()
            : type(Unknown)
        {};


        /** Explicit constructor
         *
         * @param vt   ValueType for the new TypeID. If you want to use the
         *             next available ValueType, generate it using the
         *             'newValueType()' method, as indicated in the example in
         *             the documentation.
         */
        TypeID(ValueType vt)
            : type(vt)
        {};


        /** Explicit constructor
         *
         * @param name string name for ValueType, first search tString and
         *             then search mapUserTypeID, if it is not found, create
         *             it.
         */
        TypeID(std::string name);


        /// Equality requires all fields to be the same
        virtual bool operator==(const TypeID& right) const
        { return type==right.type; };


        /// This ordering is somewhat arbitrary but is required to be able
        /// to use an TypeID as an index to a std::map. If an application
        /// needs some other ordering, inherit and override this function.
        virtual bool operator<(const TypeID& right) const
        { return type < right.type; };


        /// Inequality operator
        bool operator!=(const TypeID& right) const
        { return !(operator==(right)); };


        /// Greater than operator
        bool operator>(const TypeID& right) const
        {  return (!operator<(right) && !operator==(right)); };


        /// Less than or equal operator
        bool operator<=(const TypeID& right) const
        { return (operator<(right) || operator==(right)); };


        /// Greater than or equal operator
        bool operator>=(const TypeID& right) const
        { return !(operator<(right)); };


        /// Assignment operator
        virtual TypeID operator=(const TypeID& right);


        /// Convenience output method
        virtual std::ostream& dump(std::ostream& s) const;


        /// Return true if this is a valid TypeID. Basically just
        /// checks that the enum is defined
        virtual bool isValid() const;


        /// Destructor
        virtual ~TypeID() {};


        /** Static method to add new TypeID's
         * @param s      Identifying string for the new TypeID
         */
        static ValueType newValueType(const std::string& s);


        /// Type of the value
        ValueType type;


        /// Map holding type descriptions
        static std::map< ValueType, std::string > tStrings;


    public:

        class Initializer
        {
        public:
            Initializer();
        };

        static Initializer TypeIDsingleton;

    public:

        /** Static method to get the user registered TypeID by name string
         * @param name      Identifying string for the new TypeID
         * @return          The desired TypeID
         */
        static TypeID byName(std::string name)
            throw(InvalidRequest);

        /** Static method to add new TypeID's by name string
         * @param name      Identifying string for the new TypeID
         * @param desc      Descriptions of the new TypeID
         * @return          The new TypeID
         */
        static TypeID regByName(std::string name,std::string desc);

        /// unregister a TypeID by it's name string
        static void unregByName(std::string name);

        /// unregister all TypeIDs registered by name string
        static void unregAll();

    private:

        /// Have user deined TypeIDs been registered ?
        static bool bUserTypeIDRegistered;

        /// Map holding user defined TypeIDs by a string
        static std::map<std::string,TypeID> mapUserTypeID;


    }; // End of class 'TypeID'



    namespace StringUtils
    {
        /// convert this object to a string representation
        std::string asString(const TypeID& p);
    }



    /// stream output for TypeID
    std::ostream& operator<<(std::ostream& s, const TypeID& p);


    bool IsCarrierPhase(const RinexObsType& rot);

    int GetCarrierBand(const RinexObsType& rot);

    int GetCarrierBand(const RinexObsID& roi);

    TypeID::ValueType ConvertToTypeID(const RinexObsType& rot,
                                      const RinexSatID& sat);

    TypeID::ValueType ConvertToTypeID(const RinexObsID& roi,
                                      const RinexSatID& sat);

    TypeID::ValueType ConvertToTypeID(std::string &str);

}  // End of namespace gpstk

#endif   // GPSTK_TYPEID_HPP
