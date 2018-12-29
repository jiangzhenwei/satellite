#pragma ident "$Id$"

/**
 * @file LinearCombinations.cpp
 * This class defines handy linear combinations of GDS data.
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
//  Dagoberto Salazar - gAGE ( http://www.gage.es ). 2007, 2008, 2009
//
//============================================================================


#include "LinearCombinations.hpp"


namespace gpstk
{

    LinearCombinations::LinearCombinations()
    {
         // Coefficients of GPS Q combinations, L1 and L2
        const double xg12( ( GAMMA_GPS_L1L2 + 1.0 )/( GAMMA_GPS_L1L2 - 1.0 ) );
        const double yg12( ( 2.0                  )/( GAMMA_GPS_L1L2 - 1.0 ) );
        const double zg12( ( 2.0 * GAMMA_GPS_L1L2 )/( GAMMA_GPS_L1L2 - 1.0 ) );

         // Coefficients of GPS Q combinations, L1 and L5
        const double xg15( ( GAMMA_GPS_L1L5 + 1.0 )/( GAMMA_GPS_L1L5 - 1.0 ) );
        const double yg15( ( 2.0                  )/( GAMMA_GPS_L1L5 - 1.0 ) );
        const double zg15( ( 2.0 * GAMMA_GPS_L1L5 )/( GAMMA_GPS_L1L5 - 1.0 ) );


         // Coefficients of Galileo Q combinations, E1 and E5a
        const double xe( ( GAMMA_GAL_L1L5 + 1.0 )/( GAMMA_GAL_L1L5 - 1.0 ) );
        const double ye( ( 2.0                  )/( GAMMA_GAL_L1L5 - 1.0 ) );
        const double ze( ( 2.0 * GAMMA_GAL_L1L5 )/( GAMMA_GAL_L1L5 - 1.0 ) );


         // Coefficients of BDS Q combinations, B1 and B2
        const double xc( ( GAMMA_BDS_L2L7 + 1.0 )/( GAMMA_BDS_L2L7 - 1.0 ) );
        const double yc( ( 2.0                  )/( GAMMA_BDS_L2L7 - 1.0 ) );
        const double zc( ( 2.0 * GAMMA_BDS_L2L7 )/( GAMMA_BDS_L2L7 - 1.0 ) );


         // Definition to compute the code with minus ionospheric delays
         // in the GPS L1 frequency
        q1CombOfGPSL1L2.header              =   TypeID::Q1OfL1L2;
        q1CombOfGPSL1L2.body[TypeID::C1]    =   +xg12;
        q1CombOfGPSL1L2.body[TypeID::C2]    =   -yg12;

         // Definition to compute the code with minus ionospheric delays
         // in the GPS L2 frequency
        q2CombOfGPSL1L2.header              =   TypeID::Q2OfL1L2;
        q2CombOfGPSL1L2.body[TypeID::C1]    =   +zg12;
        q2CombOfGPSL1L2.body[TypeID::C2]    =   -xg12;


         // Definition to compute the code with minus ionospheric delays
         // in the GPS L1 frequency
        q1CombOfGPSL1L5.header              =   TypeID::Q1OfL1L5;
        q1CombOfGPSL1L5.body[TypeID::C1]    =   +xg15;
        q1CombOfGPSL1L5.body[TypeID::C5]    =   -yg15;

         // Definition to compute the code with minus ionospheric delays
         // in the GPS L5 frequency
        q5CombOfGPSL1L5.header              =   TypeID::Q5OfL1L5;
        q5CombOfGPSL1L5.body[TypeID::C1]    =   +zg15;
        q5CombOfGPSL1L5.body[TypeID::C5]    =   -xg15;


         // Definition to compute the code with minus ionospheric delays

         // Definition to compute the code with minus ionospheric delays
         // in the Galileo L1 frequency
        q1CombOfGALL1L5.header              =   TypeID::Q1OfL1L5;
        q1CombOfGALL1L5.body[TypeID::C1]    =   +xe;
        q1CombOfGALL1L5.body[TypeID::C5]    =   -ye;

         // Definition to compute the code with minus ionospheric delays
         // in the Galileo L5 frequency
        q5CombOfGALL1L5.header              =   TypeID::Q5OfL1L5;
        q5CombOfGALL1L5.body[TypeID::C1]    =   +ze;
        q5CombOfGALL1L5.body[TypeID::C5]    =   -xe;


         // Definition to compute the code with minus ionospheric delays
         // in the BDS L2 frequency
        q2CombOfBDSL2L7.header              =   TypeID::Q2OfL2L7;
        q2CombOfBDSL2L7.body[TypeID::C2]    =   +xc;
        q2CombOfBDSL2L7.body[TypeID::C7]    =   -yc;

         // Definition to compute the code with minus ionospheric delays
         // in the BDS L7 frequency
        q7CombOfBDSL2L7.header              =   TypeID::Q7OfL2L7;
        q7CombOfBDSL2L7.body[TypeID::C2]    =   +zc;
        q7CombOfBDSL2L7.body[TypeID::C7]    =   -xc;




        // Coefficients of GPS PC/LC combinations, L1 and L2
        const double ag12( GAMMA_GPS_L1L2/( GAMMA_GPS_L1L2 - 1.0 ) );
        const double bg12( 1.0           /( GAMMA_GPS_L1L2 - 1.0 ) );


        // Coefficients of GPS PC/LC combinations, L1 and L5
        const double ag15( GAMMA_GPS_L1L5/( GAMMA_GPS_L1L5 - 1.0 ) );
        const double bg15( 1.0           /( GAMMA_GPS_L1L5 - 1.0 ) );


        // Coefficients of Galileo PC/LC combinations, L1 and L5
        const double ae15( GAMMA_GAL_L1L5/( GAMMA_GAL_L1L5 - 1.0 ) );
        const double be15( 1.0           /( GAMMA_GAL_L1L5 - 1.0 ) );


        // Coefficients of BDS PC/LC combinations, L2 and L7
        const double ac27( GAMMA_BDS_L2L7/( GAMMA_BDS_L2L7 - 1.0 ) );
        const double bc27( 1.0           /( GAMMA_BDS_L2L7 - 1.0 ) );



         // Definition to compute PC combination of GPS, C1 and C2
        pc12CombOfGPS.header                =   TypeID::PC12;
        pc12CombOfGPS.body[TypeID::C1]      =   +ag12;
        pc12CombOfGPS.body[TypeID::C2]      =   -bg12;

         // Definition to compute LC combination of GPS, L1 and L2
        lc12CombOfGPS.header                =   TypeID::LC12;
        lc12CombOfGPS.body[TypeID::L1]      =   +ag12;
        lc12CombOfGPS.body[TypeID::L2]      =   -bg12;


         // Definition to compute PC combination of GPS, C1 and C5
        pc15CombOfGPS.header                =   TypeID::PC15;
        pc15CombOfGPS.body[TypeID::C1]      =   +ag15;
        pc15CombOfGPS.body[TypeID::C5]      =   -bg15;

         // Definition to compute LC combination of GPS, L1 and L5
        lc15CombOfGPS.header                =   TypeID::LC15;
        lc15CombOfGPS.body[TypeID::L1]      =   +ag15;
        lc15CombOfGPS.body[TypeID::L5]      =   -bg15;


         // Definition to compute PC combination of Galileo, C1 and C5
        pc15CombOfGAL.header                =   TypeID::PC15;
        pc15CombOfGAL.body[TypeID::C1]      =   +ae15;
        pc15CombOfGAL.body[TypeID::C5]      =   -be15;

         // Definition to compute LC combination of Galileo, L1 and L5
        lc15CombOfGAL.header                =   TypeID::LC15;
        lc15CombOfGAL.body[TypeID::L1]      =   +ae15;
        lc15CombOfGAL.body[TypeID::L5]      =   -be15;


         // Definition to compute PC combination of BDS, C2 and C7
        pc27CombOfBDS.header                =   TypeID::PC27;
        pc27CombOfBDS.body[TypeID::C2]      =   +ac27;
        pc27CombOfBDS.body[TypeID::C7]      =   -bc27;

         // Definition to compute LC combination of BDS, L2 and L7
        lc27CombOfBDS.header                =   TypeID::LC27;
        lc27CombOfBDS.body[TypeID::L2]      =   +ac27;
        lc27CombOfBDS.body[TypeID::L7]      =   -bc27;




         // Definition to compute PI combination of GPS, C1 and C2
        pi12CombOfGPS.header                =   TypeID::PI12;
        pi12CombOfGPS.body[TypeID::C1]      =   -1.0;
        pi12CombOfGPS.body[TypeID::C2]      =   +1.0;

         // Definition to compute LI combination of GPS, L1 and L2
        li12CombOfGPS.header                =   TypeID::LI12;
        li12CombOfGPS.body[TypeID::L1]      =   +1.0;
        li12CombOfGPS.body[TypeID::L2]      =   -1.0;


         // Definition to compute PI combination of GPS, C1 and C5
        pi15CombOfGPS.header                =   TypeID::PI15;
        pi15CombOfGPS.body[TypeID::C1]      =   -1.0;
        pi15CombOfGPS.body[TypeID::C5]      =   +1.0;

         // Definition to compute LI combination of GPS, L1 and L5
        li15CombOfGPS.header                =   TypeID::LI15;
        li15CombOfGPS.body[TypeID::L1]      =   +1.0;
        li15CombOfGPS.body[TypeID::L5]      =   -1.0;


         // Definition to compute PI combination of Galileo, C1 and C5
        pi15CombOfGAL.header                =   TypeID::PI15;
        pi15CombOfGAL.body[TypeID::C1]      =   -1.0;
        pi15CombOfGAL.body[TypeID::C5]      =   +1.0;

         // Definition to compute LI combination of Galileo, L1 and L5
        li15CombOfGAL.header                =   TypeID::LI15;
        li15CombOfGAL.body[TypeID::L1]      =   +1.0;
        li15CombOfGAL.body[TypeID::L5]      =   -1.0;


         // Definition to compute PI combination of BDS, C2 and C7
        pi27CombOfBDS.header                =   TypeID::PI27;
        pi27CombOfBDS.body[TypeID::C2]      =   -1.0;
        pi27CombOfBDS.body[TypeID::C7]      =   +1.0;

         // Definition to compute LI combination of BDS, L2 and L7
        li27CombOfBDS.header                =   TypeID::LI27;
        li27CombOfBDS.body[TypeID::L2]      =   +1.0;
        li27CombOfBDS.body[TypeID::L7]      =   -1.0;




         // Coefficients of GPS PW/LW/PN/LN combinations, L1 and L2
        double cg12( L1_FREQ_GPS/(L1_FREQ_GPS - L2_FREQ_GPS) );
        double dg12( L2_FREQ_GPS/(L1_FREQ_GPS - L2_FREQ_GPS) );
        double eg12( L1_FREQ_GPS/(L1_FREQ_GPS + L2_FREQ_GPS) );
        double fg12( L2_FREQ_GPS/(L1_FREQ_GPS + L2_FREQ_GPS) );

         // Coefficients of GPS PW/LW/PN/LN combinations, L1 and L5
        double cg15( L1_FREQ_GPS/(L1_FREQ_GPS - L5_FREQ_GPS) );
        double dg15( L5_FREQ_GPS/(L1_FREQ_GPS - L5_FREQ_GPS) );
        double eg15( L1_FREQ_GPS/(L1_FREQ_GPS + L5_FREQ_GPS) );
        double fg15( L5_FREQ_GPS/(L1_FREQ_GPS + L5_FREQ_GPS) );

         // Coefficients of Galileo PW/LW/PN/LN combinations, L1 and L5
        double ce15( L1_FREQ_GAL/(L1_FREQ_GAL - L5_FREQ_GAL) );
        double de15( L5_FREQ_GAL/(L1_FREQ_GAL - L5_FREQ_GAL) );
        double ee15( L1_FREQ_GAL/(L1_FREQ_GAL + L5_FREQ_GAL) );
        double fe15( L5_FREQ_GAL/(L1_FREQ_GAL + L5_FREQ_GAL) );

         // Coefficients of BDS PW/LW/PN/LN combinations, L2 and L7
        double cc27( L2_FREQ_BDS/(L2_FREQ_BDS - L7_FREQ_BDS) );
        double dc27( L7_FREQ_BDS/(L2_FREQ_BDS - L7_FREQ_BDS) );
        double ec27( L2_FREQ_BDS/(L2_FREQ_BDS + L7_FREQ_BDS) );
        double fc27( L7_FREQ_BDS/(L2_FREQ_BDS + L7_FREQ_BDS) );


         // Definition to compute PW combination of GPS, C1 and C2
        pw12CombOfGPS.header                =   TypeID::PW12;
        pw12CombOfGPS.body[TypeID::C1]      =   +cg12;
        pw12CombOfGPS.body[TypeID::C2]      =   -dg12;

         // Definition to compute LW combination of GPS, L1 and L2
        lw12CombOfGPS.header                =   TypeID::LW12;
        lw12CombOfGPS.body[TypeID::L1]      =   +cg12;
        lw12CombOfGPS.body[TypeID::L2]      =   -dg12;


         // Definition to compute PW combination of GPS, C1 and C5
        pw15CombOfGPS.header                =   TypeID::PW15;
        pw15CombOfGPS.body[TypeID::C1]      =   +cg15;
        pw15CombOfGPS.body[TypeID::C5]      =   -dg15;

         // Definition to compute LW combination of GPS, L1 and L5
        lw15CombOfGPS.header                =   TypeID::LW15;
        lw15CombOfGPS.body[TypeID::L1]      =   +cg15;
        lw15CombOfGPS.body[TypeID::L5]      =   -dg15;


         // Definition to compute PW combination of Galileo, C1 and C5
        pw15CombOfGAL.header                =   TypeID::PW15;
        pw15CombOfGAL.body[TypeID::C1]      =   +ce15;
        pw15CombOfGAL.body[TypeID::C5]      =   -de15;

         // Definition to compute LW combination of Galileo, L1 and L5
        lw15CombOfGAL.header                =   TypeID::LW15;
        lw15CombOfGAL.body[TypeID::L1]      =   +ce15;
        lw15CombOfGAL.body[TypeID::L5]      =   -de15;


         // Definition to compute PW combination of BDS, C2 and C7
        pw27CombOfBDS.header                =   TypeID::PW27;
        pw27CombOfBDS.body[TypeID::C2]      =   +cc27;
        pw27CombOfBDS.body[TypeID::C7]      =   -dc27;

         // Definition to compute LW combination of BDS, L2 and L7
        lw27CombOfBDS.header                =   TypeID::LW27;
        lw27CombOfBDS.body[TypeID::L2]      =   +cc27;
        lw27CombOfBDS.body[TypeID::L7]      =   -dc27;




         // Definition to compute PN combination of GPS, C1 and C2
        pn12CombOfGPS.header                =   TypeID::PN12;
        pn12CombOfGPS.body[TypeID::C1]      =   +eg12;
        pn12CombOfGPS.body[TypeID::C2]      =   +fg12;

         // Definition to compute LN combination of GPS, L1 and L2
        ln12CombOfGPS.header                =   TypeID::LN12;
        ln12CombOfGPS.body[TypeID::L1]      =   +eg12;
        ln12CombOfGPS.body[TypeID::L2]      =   +fg12;


         // Definition to compute PN combination of GPS, C1 and C5
        pn15CombOfGPS.header                =   TypeID::PN15;
        pn15CombOfGPS.body[TypeID::C1]      =   +eg15;
        pn15CombOfGPS.body[TypeID::C5]      =   +fg15;

         // Definition to compute LN combination of GPS, L1 and L5
        ln15CombOfGPS.header                =   TypeID::LN15;
        ln15CombOfGPS.body[TypeID::L1]      =   +eg15;
        ln15CombOfGPS.body[TypeID::L2]      =   +fg15;


         // Definition to compute PN combination of Galileo, C1 and C5
        pn15CombOfGAL.header                =   TypeID::PN15;
        pn15CombOfGAL.body[TypeID::C1]      =   +ee15;
        pn15CombOfGAL.body[TypeID::C5]      =   +fe15;

         // Definition to compute LN combination of Galileo, L1 and L5
        ln15CombOfGAL.header                =   TypeID::LN15;
        ln15CombOfGAL.body[TypeID::L1]      =   +ee15;
        ln15CombOfGAL.body[TypeID::L5]      =   +fe15;


         // Definition to compute PN combination of BDS, C2 and C7
        pn27CombOfBDS.header                =   TypeID::PN27;
        pn27CombOfBDS.body[TypeID::C2]      =   +ec27;
        pn27CombOfBDS.body[TypeID::C7]      =   +fc27;

         // Definition to compute LN combination of BDS, C2 and C7
        ln27CombOfBDS.header                =   TypeID::LN27;
        ln27CombOfBDS.body[TypeID::L2]      =   +ec27;
        ln27CombOfBDS.body[TypeID::L7]      =   +fc27;




         // Definition to compute the MW combination of GPS, L1 and L2
        mw12CombOfGPS.header                =   TypeID::MW12;
        mw12CombOfGPS.body[TypeID::L1]      =   +cg12;
        mw12CombOfGPS.body[TypeID::L2]      =   -dg12;
        mw12CombOfGPS.body[TypeID::C1]      =   -eg12;
        mw12CombOfGPS.body[TypeID::C2]      =   -fg12;

         // Definition to compute the MW combination of GPS, L1 and L5
        mw15CombOfGPS.header                =   TypeID::MW15;
        mw15CombOfGPS.body[TypeID::L1]      =   +cg15;
        mw15CombOfGPS.body[TypeID::L5]      =   -dg15;
        mw15CombOfGPS.body[TypeID::C1]      =   -eg15;
        mw15CombOfGPS.body[TypeID::C5]      =   -fg15;

         // Definition to compute the MW combination of Galileo, L1 and L5
        mw15CombOfGAL.header                =   TypeID::MW15;
        mw15CombOfGAL.body[TypeID::L1]      =   +ce15;
        mw15CombOfGAL.body[TypeID::L5]      =   -de15;
        mw15CombOfGAL.body[TypeID::C1]      =   -ee15;
        mw15CombOfGAL.body[TypeID::C5]      =   -fe15;

         // Definition to compute the MW combination of BDS, L2 and L7
        mw27CombOfBDS.header                =   TypeID::MW27;
        mw27CombOfBDS.body[TypeID::L2]      =   +cc27;
        mw27CombOfBDS.body[TypeID::L7]      =   -dc27;
        mw27CombOfBDS.body[TypeID::C2]      =   -ec27;
        mw27CombOfBDS.body[TypeID::C7]      =   -fc27;




         // Definition to compute prefit residual of GPS C1
        c1PrefitOfGPS.header                    =   TypeID::prefitC1;
        c1PrefitOfGPS.body[TypeID::C1]          =   +1.0;
        c1PrefitOfGPS.body[TypeID::rho]         =   -1.0;
        c1PrefitOfGPS.body[TypeID::cdtSat]      =   +1.0;
        c1PrefitOfGPS.body[TypeID::relativity]  =   -1.0;
        c1PrefitOfGPS.body[TypeID::gravDelay]   =   -1.0;
        c1PrefitOfGPS.body[TypeID::tropoSlant]  =   -1.0;
        c1PrefitOfGPS.body[TypeID::ionoL1]      =   -1.0;

         // Definition to compute prefit residual of GPS C2
        c2PrefitOfGPS.header                    =   TypeID::prefitC2;
        c2PrefitOfGPS.body[TypeID::C2]          =   +1.0;
        c2PrefitOfGPS.body[TypeID::rho]         =   -1.0;
        c2PrefitOfGPS.body[TypeID::cdtSat]      =   +1.0;
        c2PrefitOfGPS.body[TypeID::relativity]  =   -1.0;
        c2PrefitOfGPS.body[TypeID::gravDelay]   =   -1.0;
        c2PrefitOfGPS.body[TypeID::tropoSlant]  =   -1.0;
        c2PrefitOfGPS.body[TypeID::ionoL1]      =   -GAMMA_GPS_L1L2;

         // Definition to compute prefit residual of GPS C5
        c5PrefitOfGPS.header                    =   TypeID::prefitC5;
        c5PrefitOfGPS.body[TypeID::C2]          =   +1.0;
        c5PrefitOfGPS.body[TypeID::rho]         =   -1.0;
        c5PrefitOfGPS.body[TypeID::cdtSat]      =   +1.0;
        c5PrefitOfGPS.body[TypeID::relativity]  =   -1.0;
        c5PrefitOfGPS.body[TypeID::gravDelay]   =   -1.0;
        c5PrefitOfGPS.body[TypeID::tropoSlant]  =   -1.0;
        c5PrefitOfGPS.body[TypeID::ionoL1]      =   -GAMMA_GPS_L1L5;


         // Definition to compute prefit residual of Galileo C1
        c1PrefitOfGAL.header                    =   TypeID::prefitC1;
        c1PrefitOfGAL.body[TypeID::C1]          =   +1.0;
        c1PrefitOfGAL.body[TypeID::rho]         =   -1.0;
        c1PrefitOfGAL.body[TypeID::cdtSat]      =   +1.0;
        c1PrefitOfGAL.body[TypeID::relativity]  =   -1.0;
        c1PrefitOfGAL.body[TypeID::gravDelay]   =   -1.0;
        c1PrefitOfGAL.body[TypeID::tropoSlant]  =   -1.0;
        c1PrefitOfGAL.body[TypeID::ionoL1]      =   -1.0;

         // Definition to compute prefit residual of Galileo C5
        c5PrefitOfGAL.header                    =   TypeID::prefitC5;
        c5PrefitOfGAL.body[TypeID::C5]          =   +1.0;
        c5PrefitOfGAL.body[TypeID::rho]         =   -1.0;
        c5PrefitOfGAL.body[TypeID::cdtSat]      =   +1.0;
        c5PrefitOfGAL.body[TypeID::relativity]  =   -1.0;
        c5PrefitOfGAL.body[TypeID::gravDelay]   =   -1.0;
        c5PrefitOfGAL.body[TypeID::tropoSlant]  =   -1.0;
        c5PrefitOfGAL.body[TypeID::ionoL1]      =   -GAMMA_GAL_L1L5;

         // Definition to compute prefit residual of Galileo C6
        c6PrefitOfGAL.header                    =   TypeID::prefitC6;
        c6PrefitOfGAL.body[TypeID::C6]          =   +1.0;
        c6PrefitOfGAL.body[TypeID::rho]         =   -1.0;
        c6PrefitOfGAL.body[TypeID::cdtSat]      =   +1.0;
        c6PrefitOfGAL.body[TypeID::relativity]  =   -1.0;
        c6PrefitOfGAL.body[TypeID::gravDelay]   =   -1.0;
        c6PrefitOfGAL.body[TypeID::tropoSlant]  =   -1.0;
        c6PrefitOfGAL.body[TypeID::ionoL1]      =   -GAMMA_GAL_L1L6;

         // Definition to compute prefit residual of Galileo C7
        c7PrefitOfGAL.header                    =   TypeID::prefitC7;
        c7PrefitOfGAL.body[TypeID::C7]          =   +1.0;
        c7PrefitOfGAL.body[TypeID::rho]         =   -1.0;
        c7PrefitOfGAL.body[TypeID::cdtSat]      =   +1.0;
        c7PrefitOfGAL.body[TypeID::relativity]  =   -1.0;
        c7PrefitOfGAL.body[TypeID::gravDelay]   =   -1.0;
        c7PrefitOfGAL.body[TypeID::tropoSlant]  =   -1.0;
        c7PrefitOfGAL.body[TypeID::ionoL1]      =   -GAMMA_GAL_L1L7;

         // Definition to compute prefit residual of Galileo C6
        c8PrefitOfGAL.header                    =   TypeID::prefitC8;
        c8PrefitOfGAL.body[TypeID::C8]          =   +1.0;
        c8PrefitOfGAL.body[TypeID::rho]         =   -1.0;
        c8PrefitOfGAL.body[TypeID::cdtSat]      =   +1.0;
        c8PrefitOfGAL.body[TypeID::relativity]  =   -1.0;
        c8PrefitOfGAL.body[TypeID::gravDelay]   =   -1.0;
        c8PrefitOfGAL.body[TypeID::tropoSlant]  =   -1.0;
        c8PrefitOfGAL.body[TypeID::ionoL1]      =   -GAMMA_GAL_L1L8;


         // Definition to compute prefit residual of BDS C2
        c2PrefitOfBDS.header                    =   TypeID::prefitC2;
        c2PrefitOfBDS.body[TypeID::C2]          =   +1.0;
        c2PrefitOfBDS.body[TypeID::rho]         =   -1.0;
        c2PrefitOfBDS.body[TypeID::cdtSat]      =   +1.0;
        c2PrefitOfBDS.body[TypeID::relativity]  =   -1.0;
        c2PrefitOfBDS.body[TypeID::gravDelay]   =   -1.0;
        c2PrefitOfBDS.body[TypeID::tropoSlant]  =   -1.0;
        c2PrefitOfBDS.body[TypeID::ionoL2]      =   -1.0;

         // Definition to compute prefit residual of BDS C7
        c7PrefitOfBDS.header                    =   TypeID::prefitC7;
        c7PrefitOfBDS.body[TypeID::C7]          =   +1.0;
        c7PrefitOfBDS.body[TypeID::rho]         =   -1.0;
        c7PrefitOfBDS.body[TypeID::cdtSat]      =   +1.0;
        c7PrefitOfBDS.body[TypeID::relativity]  =   -1.0;
        c7PrefitOfBDS.body[TypeID::gravDelay]   =   -1.0;
        c7PrefitOfBDS.body[TypeID::tropoSlant]  =   -1.0;
        c7PrefitOfBDS.body[TypeID::ionoL2]      =   -GAMMA_BDS_L2L7;

         // Definition to compute prefit residual of Galileo C6
        c6PrefitOfBDS.header                    =   TypeID::prefitC6;
        c6PrefitOfBDS.body[TypeID::C6]          =   +1.0;
        c6PrefitOfBDS.body[TypeID::rho]         =   -1.0;
        c6PrefitOfBDS.body[TypeID::cdtSat]      =   +1.0;
        c6PrefitOfBDS.body[TypeID::relativity]  =   -1.0;
        c6PrefitOfBDS.body[TypeID::gravDelay]   =   -1.0;
        c6PrefitOfBDS.body[TypeID::tropoSlant]  =   -1.0;
        c6PrefitOfBDS.body[TypeID::ionoL2]      =   -GAMMA_BDS_L2L6;



         // Definition to compute prefit residual of GPS L1
        l1PrefitOfGPS.header                      =   TypeID::prefitL1;
        l1PrefitOfGPS.body[TypeID::L1]            =   +1.0;
        l1PrefitOfGPS.body[TypeID::rho]           =   -1.0;
        l1PrefitOfGPS.body[TypeID::cdtSat]        =   +1.0;
        l1PrefitOfGPS.body[TypeID::relativity]    =   -1.0;
        l1PrefitOfGPS.body[TypeID::gravDelay]     =   -1.0;
        l1PrefitOfGPS.body[TypeID::tropoSlant]    =   -1.0;
        l1PrefitOfGPS.body[TypeID::ionoL1]        =   +1.0;
        l1PrefitOfGPS.body[TypeID::windUp]        =   -L1_WAVELENGTH_GPS/TWO_PI;

         // Definition to compute prefit residual of GPS L2
        l2PrefitOfGPS.header                      =   TypeID::prefitL2;
        l2PrefitOfGPS.body[TypeID::L2]            =   +1.0;
        l2PrefitOfGPS.body[TypeID::rho]           =   -1.0;
        l2PrefitOfGPS.body[TypeID::cdtSat]        =   +1.0;
        l2PrefitOfGPS.body[TypeID::relativity]    =   -1.0;
        l2PrefitOfGPS.body[TypeID::gravDelay]     =   -1.0;
        l2PrefitOfGPS.body[TypeID::tropoSlant]    =   -1.0;
        l2PrefitOfGPS.body[TypeID::ionoL1]        =   +GAMMA_GPS_L1L2;
        l2PrefitOfGPS.body[TypeID::windUp]        =   -L2_WAVELENGTH_GPS/TWO_PI;

         // Definition to compute prefit residual of GPS L5
        l5PrefitOfGPS.header                      =   TypeID::prefitL5;
        l5PrefitOfGPS.body[TypeID::L5]            =   +1.0;
        l5PrefitOfGPS.body[TypeID::rho]           =   -1.0;
        l5PrefitOfGPS.body[TypeID::cdtSat]        =   +1.0;
        l5PrefitOfGPS.body[TypeID::relativity]    =   -1.0;
        l5PrefitOfGPS.body[TypeID::gravDelay]     =   -1.0;
        l5PrefitOfGPS.body[TypeID::tropoSlant]    =   -1.0;
        l5PrefitOfGPS.body[TypeID::ionoL1]        =   +GAMMA_GPS_L1L5;
        l5PrefitOfGPS.body[TypeID::windUp]        =   -L5_WAVELENGTH_GPS/TWO_PI;


         // Definition to compute prefit residual of Galileo L1
        l1PrefitOfGAL.header                      =   TypeID::prefitL1;
        l1PrefitOfGAL.body[TypeID::L1]            =   +1.0;
        l1PrefitOfGAL.body[TypeID::rho]           =   -1.0;
        l1PrefitOfGAL.body[TypeID::cdtSat]        =   +1.0;
        l1PrefitOfGAL.body[TypeID::relativity]    =   -1.0;
        l1PrefitOfGAL.body[TypeID::gravDelay]     =   -1.0;
        l1PrefitOfGAL.body[TypeID::tropoSlant]    =   -1.0;
        l1PrefitOfGAL.body[TypeID::ionoL1]        =   +1.0;
        l1PrefitOfGAL.body[TypeID::windUp]        =   -L1_WAVELENGTH_GAL/TWO_PI;

         // Definition to compute prefit residual of Galileo L5
        l5PrefitOfGAL.header                      =   TypeID::prefitL5;
        l5PrefitOfGAL.body[TypeID::L5]            =   +1.0;
        l5PrefitOfGAL.body[TypeID::rho]           =   -1.0;
        l5PrefitOfGAL.body[TypeID::cdtSat]        =   +1.0;
        l5PrefitOfGAL.body[TypeID::relativity]    =   -1.0;
        l5PrefitOfGAL.body[TypeID::gravDelay]     =   -1.0;
        l5PrefitOfGAL.body[TypeID::tropoSlant]    =   -1.0;
        l5PrefitOfGAL.body[TypeID::ionoL5]        =   +GAMMA_GAL_L1L5;
        l5PrefitOfGAL.body[TypeID::windUp]        =   -L5_WAVELENGTH_GAL/TWO_PI;

         // Definition to compute prefit residual of Galileo L6
        l6PrefitOfGAL.header                      =   TypeID::prefitL6;
        l6PrefitOfGAL.body[TypeID::L6]            =   +1.0;
        l6PrefitOfGAL.body[TypeID::rho]           =   -1.0;
        l6PrefitOfGAL.body[TypeID::cdtSat]        =   +1.0;
        l6PrefitOfGAL.body[TypeID::relativity]    =   -1.0;
        l6PrefitOfGAL.body[TypeID::gravDelay]     =   -1.0;
        l6PrefitOfGAL.body[TypeID::tropoSlant]    =   -1.0;
        l6PrefitOfGAL.body[TypeID::ionoL1]        =   +GAMMA_GAL_L1L6;
        l6PrefitOfGAL.body[TypeID::windUp]        =   -L6_WAVELENGTH_GAL/TWO_PI;

         // Definition to compute prefit residual of Galileo L7
        l7PrefitOfGAL.header                      =   TypeID::prefitL7;
        l7PrefitOfGAL.body[TypeID::L7]            =   +1.0;
        l7PrefitOfGAL.body[TypeID::rho]           =   -1.0;
        l7PrefitOfGAL.body[TypeID::cdtSat]        =   +1.0;
        l7PrefitOfGAL.body[TypeID::relativity]    =   -1.0;
        l7PrefitOfGAL.body[TypeID::gravDelay]     =   -1.0;
        l7PrefitOfGAL.body[TypeID::tropoSlant]    =   -1.0;
        l7PrefitOfGAL.body[TypeID::ionoL1]        =   +GAMMA_GAL_L1L7;
        l7PrefitOfGAL.body[TypeID::windUp]        =   -L7_WAVELENGTH_GAL/TWO_PI;

         // Definition to compute prefit residual of Galileo L8
        l8PrefitOfGAL.header                      =   TypeID::prefitL8;
        l8PrefitOfGAL.body[TypeID::L8]            =   +1.0;
        l8PrefitOfGAL.body[TypeID::rho]           =   -1.0;
        l8PrefitOfGAL.body[TypeID::cdtSat]        =   +1.0;
        l8PrefitOfGAL.body[TypeID::relativity]    =   -1.0;
        l8PrefitOfGAL.body[TypeID::gravDelay]     =   -1.0;
        l8PrefitOfGAL.body[TypeID::tropoSlant]    =   -1.0;
        l8PrefitOfGAL.body[TypeID::ionoL1]        =   +GAMMA_GAL_L1L8;
        l8PrefitOfGAL.body[TypeID::windUp]        =   -L8_WAVELENGTH_GAL/TWO_PI;



         // Definition to compute prefit residual of BDS L2
        l2PrefitOfBDS.header                      =   TypeID::prefitL2;
        l2PrefitOfBDS.body[TypeID::L2]            =   +1.0;
        l2PrefitOfBDS.body[TypeID::rho]           =   -1.0;
        l2PrefitOfBDS.body[TypeID::cdtSat]        =   +1.0;
        l2PrefitOfBDS.body[TypeID::relativity]    =   -1.0;
        l2PrefitOfBDS.body[TypeID::gravDelay]     =   -1.0;
        l2PrefitOfBDS.body[TypeID::tropoSlant]    =   -1.0;
        l2PrefitOfBDS.body[TypeID::ionoL2]        =   +1.0;
        l2PrefitOfBDS.body[TypeID::windUp]        =   -L2_WAVELENGTH_BDS/TWO_PI;

         // Definition to compute prefit residual of BDS L7
        l7PrefitOfBDS.header                      =   TypeID::prefitL7;
        l7PrefitOfBDS.body[TypeID::L7]            =   +1.0;
        l7PrefitOfBDS.body[TypeID::rho]           =   -1.0;
        l7PrefitOfBDS.body[TypeID::cdtSat]        =   +1.0;
        l7PrefitOfBDS.body[TypeID::relativity]    =   -1.0;
        l7PrefitOfBDS.body[TypeID::gravDelay]     =   -1.0;
        l7PrefitOfBDS.body[TypeID::tropoSlant]    =   -1.0;
        l7PrefitOfBDS.body[TypeID::ionoL2]        =   +GAMMA_BDS_L2L7;
        l7PrefitOfBDS.body[TypeID::windUp]        =   -L7_WAVELENGTH_BDS/TWO_PI;

         // Definition to compute prefit residual of BDS L6
        l6PrefitOfBDS.header                      =   TypeID::prefitL6;
        l6PrefitOfBDS.body[TypeID::L7]            =   +1.0;
        l6PrefitOfBDS.body[TypeID::rho]           =   -1.0;
        l6PrefitOfBDS.body[TypeID::cdtSat]        =   +1.0;
        l6PrefitOfBDS.body[TypeID::relativity]    =   -1.0;
        l6PrefitOfBDS.body[TypeID::gravDelay]     =   -1.0;
        l6PrefitOfBDS.body[TypeID::tropoSlant]    =   -1.0;
        l6PrefitOfBDS.body[TypeID::ionoL2]        =   +GAMMA_BDS_L2L6;
        l6PrefitOfBDS.body[TypeID::windUp]        =   -L6_WAVELENGTH_BDS/TWO_PI;




         // Definition to compute prefit residual of GPS PC, L1 + L2
        pc12PrefitOfGPS.header                                  =   TypeID::prefitC12;
        pc12PrefitOfGPS.body[TypeID::PC12]                      =   +1.0;
        pc12PrefitOfGPS.body[TypeID::rho]                       =   -1.0;
        pc12PrefitOfGPS.body[TypeID::cdtSat]                    =   +1.0;
        pc12PrefitOfGPS.body[TypeID::relativity]                =   -1.0;
        pc12PrefitOfGPS.body[TypeID::gravDelay]                 =   -1.0;
        pc12PrefitOfGPS.body[TypeID::tropoSlant]                =   -1.0;

         // Definition to compute prefit residual of GPS PC, L1 + L2, without clock
        pc12PrefitOfGPSWithoutClock.header                      =   TypeID::prefitC12WithoutClock;
        pc12PrefitOfGPSWithoutClock.body[TypeID::PC12]          =   +1.0;
        pc12PrefitOfGPSWithoutClock.body[TypeID::rho]           =   -1.0;
        pc12PrefitOfGPSWithoutClock.body[TypeID::relativity]    =   -1.0;
        pc12PrefitOfGPSWithoutClock.body[TypeID::gravDelay]     =   -1.0;
        pc12PrefitOfGPSWithoutClock.body[TypeID::tropoSlant]    =   -1.0;

         // Definition to compute prefit residual of GPS PC, L1 + L2, with sat clock
        pc12PrefitOfGPSWithSatClock.header                      =   TypeID::prefitC12WithSatClock;
        pc12PrefitOfGPSWithSatClock.body[TypeID::PC12]          =   +1.0;
        pc12PrefitOfGPSWithSatClock.body[TypeID::rho]           =   -1.0;
        pc12PrefitOfGPSWithSatClock.body[TypeID::cdtSat]        =   +1.0;
        pc12PrefitOfGPSWithSatClock.body[TypeID::relativity]    =   -1.0;
        pc12PrefitOfGPSWithSatClock.body[TypeID::gravDelay]     =   -1.0;
        pc12PrefitOfGPSWithSatClock.body[TypeID::tropoSlant]    =   -1.0;

         // Definition to compute prefit residual of GPS PC, L1 + L2, with sta clock
        pc12PrefitOfGPSWithStaClock.header                      =   TypeID::prefitC12WithStaClock;
        pc12PrefitOfGPSWithStaClock.body[TypeID::PC12]          =   +1.0;
        pc12PrefitOfGPSWithStaClock.body[TypeID::rho]           =   -1.0;
        pc12PrefitOfGPSWithStaClock.body[TypeID::cdtSta]        =   -1.0;
        pc12PrefitOfGPSWithStaClock.body[TypeID::relativity]    =   -1.0;
        pc12PrefitOfGPSWithStaClock.body[TypeID::gravDelay]     =   -1.0;
        pc12PrefitOfGPSWithStaClock.body[TypeID::tropoSlant]    =   -1.0;

        // Definition to compute prefit residual of GPS PC, L1 + L2, for PCE
        pc12PrefitOfGPSForPCE.header                            =   TypeID::prefitC12ForPCE;
        pc12PrefitOfGPSForPCE.body[TypeID::PC12]                =   +1.0;
        pc12PrefitOfGPSForPCE.body[TypeID::rho]                 =   -1.0;
        pc12PrefitOfGPSForPCE.body[TypeID::relativity]          =   -1.0;
        pc12PrefitOfGPSForPCE.body[TypeID::gravDelay]           =   -1.0;
        pc12PrefitOfGPSForPCE.body[TypeID::tropoSlant]          =   -1.0;

        // Definition to compute prefit residual of GPS PC, L1 + L2, for POD
        pc12PrefitOfGPSForPOD.header                            =   TypeID::prefitC12ForPOD;
        pc12PrefitOfGPSForPOD.body[TypeID::PC12]                =   +1.0;
        pc12PrefitOfGPSForPOD.body[TypeID::relativity]          =   -1.0;
        pc12PrefitOfGPSForPOD.body[TypeID::gravDelay]           =   -1.0;
        pc12PrefitOfGPSForPOD.body[TypeID::tropoSlant]          =   -1.0;

        // Definition to compute prefit residual of GPS PC, L1 + L2, for PPP
        pc12PrefitOfGPSForPPP.header                            =   TypeID::prefitC12ForPPP;
        pc12PrefitOfGPSForPPP.body[TypeID::PC12]                =   +1.0;
        pc12PrefitOfGPSForPPP.body[TypeID::rho]                 =   -1.0;
        pc12PrefitOfGPSForPPP.body[TypeID::cdtSat]              =   +1.0;
        pc12PrefitOfGPSForPPP.body[TypeID::relativity]          =   -1.0;
        pc12PrefitOfGPSForPPP.body[TypeID::gravDelay]           =   -1.0;
        pc12PrefitOfGPSForPPP.body[TypeID::tropoSlant]          =   -1.0;


         // Definition to compute prefit residual of GPS LC, L1 + L2
        lc12PrefitOfGPS.header                                  =   TypeID::prefitL12;
        lc12PrefitOfGPS.body[TypeID::LC12]                      =   +1.0;
        lc12PrefitOfGPS.body[TypeID::rho]                       =   -1.0;
        lc12PrefitOfGPS.body[TypeID::cdtSat]                    =   +1.0;
        lc12PrefitOfGPS.body[TypeID::relativity]                =   -1.0;
        lc12PrefitOfGPS.body[TypeID::gravDelay]                 =   -1.0;
        lc12PrefitOfGPS.body[TypeID::tropoSlant]                =   -1.0;
        lc12PrefitOfGPS.body[TypeID::windUp]                    =   -LC_WAVELENGTH_GPS_L1L2/TWO_PI;

         // Definition to compute prefit residual of GPS LC, L1 + L2, without clock
        lc12PrefitOfGPSWithoutClock.header                      =   TypeID::prefitL12WithoutClock;
        lc12PrefitOfGPSWithoutClock.body[TypeID::LC12]          =   +1.0;
        lc12PrefitOfGPSWithoutClock.body[TypeID::rho]           =   -1.0;
        lc12PrefitOfGPSWithoutClock.body[TypeID::relativity]    =   -1.0;
        lc12PrefitOfGPSWithoutClock.body[TypeID::gravDelay]     =   -1.0;
        lc12PrefitOfGPSWithoutClock.body[TypeID::tropoSlant]    =   -1.0;
        lc12PrefitOfGPSWithoutClock.body[TypeID::windUp]        =   -LC_WAVELENGTH_GPS_L1L2/TWO_PI;

         // Definition to compute prefit residual of GPS LC, L1 + L2, with sat clock
        lc12PrefitOfGPSWithSatClock.header                      =   TypeID::prefitL12WithSatClock;
        lc12PrefitOfGPSWithSatClock.body[TypeID::LC12]          =   +1.0;
        lc12PrefitOfGPSWithSatClock.body[TypeID::rho]           =   -1.0;
        lc12PrefitOfGPSWithSatClock.body[TypeID::cdtSat]        =   +1.0;
        lc12PrefitOfGPSWithSatClock.body[TypeID::relativity]    =   -1.0;
        lc12PrefitOfGPSWithSatClock.body[TypeID::gravDelay]     =   -1.0;
        lc12PrefitOfGPSWithSatClock.body[TypeID::tropoSlant]    =   -1.0;
        lc12PrefitOfGPSWithSatClock.body[TypeID::windUp]        =   -LC_WAVELENGTH_GPS_L1L2/TWO_PI;

         // Definition to compute prefit residual of GPS LC, L1 + L2, with sta clock
        lc12PrefitOfGPSWithStaClock.header                      =   TypeID::prefitL12WithStaClock;
        lc12PrefitOfGPSWithStaClock.body[TypeID::LC12]          =   +1.0;
        lc12PrefitOfGPSWithStaClock.body[TypeID::rho]           =   -1.0;
        lc12PrefitOfGPSWithStaClock.body[TypeID::cdtSta]        =   -1.0;
        lc12PrefitOfGPSWithStaClock.body[TypeID::relativity]    =   -1.0;
        lc12PrefitOfGPSWithStaClock.body[TypeID::gravDelay]     =   -1.0;
        lc12PrefitOfGPSWithStaClock.body[TypeID::tropoSlant]    =   -1.0;
        lc12PrefitOfGPSWithStaClock.body[TypeID::windUp]        =   -LC_WAVELENGTH_GPS_L1L2/TWO_PI;

        // Definition to compute prefit residual of GPS LC, L1 + L2, for PCE
        lc12PrefitOfGPSForPCE.header                            =   TypeID::prefitL12ForPCE;
        lc12PrefitOfGPSForPCE.body[TypeID::LC12]                =   +1.0;
        lc12PrefitOfGPSForPCE.body[TypeID::rho]                 =   -1.0;
        lc12PrefitOfGPSForPCE.body[TypeID::relativity]          =   -1.0;
        lc12PrefitOfGPSForPCE.body[TypeID::gravDelay]           =   -1.0;
        lc12PrefitOfGPSForPCE.body[TypeID::tropoSlant]          =   -1.0;
        lc12PrefitOfGPSForPCE.body[TypeID::windUp]              =   -LC_WAVELENGTH_GPS_L1L2/TWO_PI;

        // Definition to compute prefit residual of GPS LC, L1 + L2, for POD
        lc12PrefitOfGPSForPOD.header                            =   TypeID::prefitL12ForPOD;
        lc12PrefitOfGPSForPOD.body[TypeID::LC12]                =   +1.0;
        lc12PrefitOfGPSForPOD.body[TypeID::relativity]          =   -1.0;
        lc12PrefitOfGPSForPOD.body[TypeID::gravDelay]           =   -1.0;
        lc12PrefitOfGPSForPOD.body[TypeID::tropoSlant]          =   -1.0;
        lc12PrefitOfGPSForPOD.body[TypeID::windUp]              =   -LC_WAVELENGTH_GPS_L1L2/TWO_PI;

        // Definition to compute prefit residual of GPS LC, L1 + L2, for PPP
        lc12PrefitOfGPSForPPP.header                            =   TypeID::prefitL12ForPPP;
        lc12PrefitOfGPSForPPP.body[TypeID::LC12]                =   +1.0;
        lc12PrefitOfGPSForPPP.body[TypeID::rho]                 =   -1.0;
        lc12PrefitOfGPSForPPP.body[TypeID::cdtSat]              =   +1.0;
        lc12PrefitOfGPSForPPP.body[TypeID::relativity]          =   -1.0;
        lc12PrefitOfGPSForPPP.body[TypeID::gravDelay]           =   -1.0;
        lc12PrefitOfGPSForPPP.body[TypeID::tropoSlant]          =   -1.0;
        lc12PrefitOfGPSForPPP.body[TypeID::windUp]              =   -LC_WAVELENGTH_GPS_L1L2/TWO_PI;


         // Definition to compute prefit residual of GPS PC, L1 + L5
        pc15PrefitOfGPS.header                                  =   TypeID::prefitC15;
        pc15PrefitOfGPS.body[TypeID::PC15]                      =   +1.0;
        pc15PrefitOfGPS.body[TypeID::rho]                       =   -1.0;
        pc15PrefitOfGPS.body[TypeID::cdtSat]                    =   +1.0;
        pc15PrefitOfGPS.body[TypeID::relativity]                =   -1.0;
        pc15PrefitOfGPS.body[TypeID::gravDelay]                 =   -1.0;
        pc15PrefitOfGPS.body[TypeID::tropoSlant]                =   -1.0;

         // Definition to compute prefit residual of GPS PC, L1 + L5, without clock
        pc15PrefitOfGPSWithoutClock.header                      =   TypeID::prefitC15WithoutClock;
        pc15PrefitOfGPSWithoutClock.body[TypeID::PC15]          =   +1.0;
        pc15PrefitOfGPSWithoutClock.body[TypeID::rho]           =   -1.0;
        pc15PrefitOfGPSWithoutClock.body[TypeID::relativity]    =   -1.0;
        pc15PrefitOfGPSWithoutClock.body[TypeID::gravDelay]     =   -1.0;
        pc15PrefitOfGPSWithoutClock.body[TypeID::tropoSlant]    =   -1.0;

         // Definition to compute prefit residual of GPS PC, L1 + L5, with sat clock
        pc15PrefitOfGPSWithSatClock.header                      =   TypeID::prefitC15WithSatClock;
        pc15PrefitOfGPSWithSatClock.body[TypeID::PC15]          =   +1.0;
        pc15PrefitOfGPSWithSatClock.body[TypeID::rho]           =   -1.0;
        pc15PrefitOfGPSWithSatClock.body[TypeID::cdtSat]        =   +1.0;
        pc15PrefitOfGPSWithSatClock.body[TypeID::relativity]    =   -1.0;
        pc15PrefitOfGPSWithSatClock.body[TypeID::gravDelay]     =   -1.0;
        pc15PrefitOfGPSWithSatClock.body[TypeID::tropoSlant]    =   -1.0;

         // Definition to compute prefit residual of GPS PC, L1 + L5, with sta clock
        pc15PrefitOfGPSWithStaClock.header                      =   TypeID::prefitC15WithStaClock;
        pc15PrefitOfGPSWithStaClock.body[TypeID::PC15]          =   +1.0;
        pc15PrefitOfGPSWithStaClock.body[TypeID::rho]           =   -1.0;
        pc15PrefitOfGPSWithStaClock.body[TypeID::cdtSta]        =   -1.0;
        pc15PrefitOfGPSWithStaClock.body[TypeID::relativity]    =   -1.0;
        pc15PrefitOfGPSWithStaClock.body[TypeID::gravDelay]     =   -1.0;
        pc15PrefitOfGPSWithStaClock.body[TypeID::tropoSlant]    =   -1.0;

        // Definition to compute prefit residual of GPS PC, L1 + L5, for PCE
        pc15PrefitOfGPSForPCE.header                            =   TypeID::prefitC15ForPCE;
        pc15PrefitOfGPSForPCE.body[TypeID::PC15]                =   +1.0;
        pc15PrefitOfGPSForPCE.body[TypeID::rho]                 =   -1.0;
        pc15PrefitOfGPSForPCE.body[TypeID::relativity]          =   -1.0;
        pc15PrefitOfGPSForPCE.body[TypeID::gravDelay]           =   -1.0;
        pc15PrefitOfGPSForPCE.body[TypeID::tropoSlant]          =   -1.0;

        // Definition to compute prefit residual of GPS PC, L1 + L5, for POD
        pc15PrefitOfGPSForPOD.header                            =   TypeID::prefitC15ForPOD;
        pc15PrefitOfGPSForPOD.body[TypeID::PC15]                =   +1.0;
        pc15PrefitOfGPSForPOD.body[TypeID::relativity]          =   -1.0;
        pc15PrefitOfGPSForPOD.body[TypeID::gravDelay]           =   -1.0;
        pc15PrefitOfGPSForPOD.body[TypeID::tropoSlant]          =   -1.0;

        // Definition to compute prefit residual of GPS PC, L1 + L5, for PPP
        pc15PrefitOfGPSForPPP.header                            =   TypeID::prefitC15ForPPP;
        pc15PrefitOfGPSForPPP.body[TypeID::PC15]                =   +1.0;
        pc15PrefitOfGPSForPPP.body[TypeID::rho]                 =   -1.0;
        pc15PrefitOfGPSForPPP.body[TypeID::cdtSat]              =   +1.0;
        pc15PrefitOfGPSForPPP.body[TypeID::relativity]          =   -1.0;
        pc15PrefitOfGPSForPPP.body[TypeID::gravDelay]           =   -1.0;
        pc15PrefitOfGPSForPPP.body[TypeID::tropoSlant]          =   -1.0;


         // Definition to compute prefit residual of GPS LC, L1 + L5
        lc15PrefitOfGPS.header                                  =   TypeID::prefitL15;
        lc15PrefitOfGPS.body[TypeID::LC15]                      =   +1.0;
        lc15PrefitOfGPS.body[TypeID::rho]                       =   -1.0;
        lc15PrefitOfGPS.body[TypeID::cdtSat]                    =   +1.0;
        lc15PrefitOfGPS.body[TypeID::relativity]                =   -1.0;
        lc15PrefitOfGPS.body[TypeID::gravDelay]                 =   -1.0;
        lc15PrefitOfGPS.body[TypeID::tropoSlant]                =   -1.0;
        lc15PrefitOfGPS.body[TypeID::windUp]                    =   -LC_WAVELENGTH_GPS_L1L5/TWO_PI;

         // Definition to compute prefit residual of GPS LC, L1 + L5, without clock
        lc15PrefitOfGPSWithoutClock.header                      =   TypeID::prefitL15WithoutClock;
        lc15PrefitOfGPSWithoutClock.body[TypeID::LC15]          =   +1.0;
        lc15PrefitOfGPSWithoutClock.body[TypeID::rho]           =   -1.0;
        lc15PrefitOfGPSWithoutClock.body[TypeID::relativity]    =   -1.0;
        lc15PrefitOfGPSWithoutClock.body[TypeID::gravDelay]     =   -1.0;
        lc15PrefitOfGPSWithoutClock.body[TypeID::tropoSlant]    =   -1.0;
        lc15PrefitOfGPSWithoutClock.body[TypeID::windUp]        =   -LC_WAVELENGTH_GPS_L1L5/TWO_PI;

         // Definition to compute prefit residual of GPS LC, L1 + L5, with sat clock
        lc15PrefitOfGPSWithSatClock.header                      =   TypeID::prefitL15WithSatClock;
        lc15PrefitOfGPSWithSatClock.body[TypeID::LC15]          =   +1.0;
        lc15PrefitOfGPSWithSatClock.body[TypeID::rho]           =   -1.0;
        lc15PrefitOfGPSWithSatClock.body[TypeID::cdtSat]        =   +1.0;
        lc15PrefitOfGPSWithSatClock.body[TypeID::relativity]    =   -1.0;
        lc15PrefitOfGPSWithSatClock.body[TypeID::gravDelay]     =   -1.0;
        lc15PrefitOfGPSWithSatClock.body[TypeID::tropoSlant]    =   -1.0;
        lc15PrefitOfGPSWithSatClock.body[TypeID::windUp]        =   -LC_WAVELENGTH_GPS_L1L5/TWO_PI;

         // Definition to compute prefit residual of GPS LC, L1 + L5, with sta clock
        lc15PrefitOfGPSWithStaClock.header                      =   TypeID::prefitL15WithStaClock;
        lc15PrefitOfGPSWithStaClock.body[TypeID::LC15]          =   +1.0;
        lc15PrefitOfGPSWithStaClock.body[TypeID::rho]           =   -1.0;
        lc15PrefitOfGPSWithStaClock.body[TypeID::cdtSta]        =   -1.0;
        lc15PrefitOfGPSWithStaClock.body[TypeID::relativity]    =   -1.0;
        lc15PrefitOfGPSWithStaClock.body[TypeID::gravDelay]     =   -1.0;
        lc15PrefitOfGPSWithStaClock.body[TypeID::tropoSlant]    =   -1.0;
        lc15PrefitOfGPSWithStaClock.body[TypeID::windUp]        =   -LC_WAVELENGTH_GPS_L1L5/TWO_PI;

        // Definition to compute prefit residual of GPS LC, L1 + L5, for PCE
        lc15PrefitOfGPSForPCE.header                            =   TypeID::prefitL15ForPCE;
        lc15PrefitOfGPSForPCE.body[TypeID::LC15]                =   +1.0;
        lc15PrefitOfGPSForPCE.body[TypeID::rho]                 =   -1.0;
        lc15PrefitOfGPSForPCE.body[TypeID::relativity]          =   -1.0;
        lc15PrefitOfGPSForPCE.body[TypeID::gravDelay]           =   -1.0;
        lc15PrefitOfGPSForPCE.body[TypeID::tropoSlant]          =   -1.0;
        lc15PrefitOfGPSForPCE.body[TypeID::windUp]              =   -LC_WAVELENGTH_GPS_L1L5/TWO_PI;

        // Definition to compute prefit residual of GPS LC, L1 + L5, for POD
        lc15PrefitOfGPSForPOD.header                            =   TypeID::prefitL15ForPOD;
        lc15PrefitOfGPSForPOD.body[TypeID::LC15]                =   +1.0;
        lc15PrefitOfGPSForPOD.body[TypeID::relativity]          =   -1.0;
        lc15PrefitOfGPSForPOD.body[TypeID::gravDelay]           =   -1.0;
        lc15PrefitOfGPSForPOD.body[TypeID::tropoSlant]          =   -1.0;
        lc15PrefitOfGPSForPOD.body[TypeID::windUp]              =   -LC_WAVELENGTH_GPS_L1L5/TWO_PI;

        // Definition to compute prefit residual of GPS LC, L1 + L5, for PPP
        lc15PrefitOfGPSForPPP.header                            =   TypeID::prefitL15ForPPP;
        lc15PrefitOfGPSForPPP.body[TypeID::LC15]                =   +1.0;
        lc15PrefitOfGPSForPPP.body[TypeID::rho]                 =   -1.0;
        lc15PrefitOfGPSForPPP.body[TypeID::cdtSat]              =   +1.0;
        lc15PrefitOfGPSForPPP.body[TypeID::relativity]          =   -1.0;
        lc15PrefitOfGPSForPPP.body[TypeID::gravDelay]           =   -1.0;
        lc15PrefitOfGPSForPPP.body[TypeID::tropoSlant]          =   -1.0;
        lc15PrefitOfGPSForPPP.body[TypeID::windUp]              =   -LC_WAVELENGTH_GPS_L1L5/TWO_PI;




         // Definition to compute prefit residual of Galileo LC, L1 + L5
        lc15PrefitOfGAL.header                                  =   TypeID::prefitL15;
        lc15PrefitOfGAL.body[TypeID::LC15]                      =   +1.0;
        lc15PrefitOfGAL.body[TypeID::rho]                       =   -1.0;
        lc15PrefitOfGAL.body[TypeID::cdtSat]                    =   +1.0;
        lc15PrefitOfGAL.body[TypeID::relativity]                =   -1.0;
        lc15PrefitOfGAL.body[TypeID::gravDelay]                 =   -1.0;
        lc15PrefitOfGAL.body[TypeID::tropoSlant]                =   -1.0;
        lc15PrefitOfGAL.body[TypeID::windUp]                    =   -LC_WAVELENGTH_GAL_L1L5/TWO_PI;

         // Definition to compute prefit residual of Galileo LC, L1 + L5, without clock
        lc15PrefitOfGALWithoutClock.header                      =   TypeID::prefitL15WithoutClock;
        lc15PrefitOfGALWithoutClock.body[TypeID::LC15]          =   +1.0;
        lc15PrefitOfGALWithoutClock.body[TypeID::rho]           =   -1.0;
        lc15PrefitOfGALWithoutClock.body[TypeID::relativity]    =   -1.0;
        lc15PrefitOfGALWithoutClock.body[TypeID::gravDelay]     =   -1.0;
        lc15PrefitOfGALWithoutClock.body[TypeID::tropoSlant]    =   -1.0;
        lc15PrefitOfGALWithoutClock.body[TypeID::windUp]        =   -LC_WAVELENGTH_GAL_L1L5/TWO_PI;

         // Definition to compute prefit residual of Galileo LC, L1 + L5, with sat clock
        lc15PrefitOfGALWithSatClock.header                      =   TypeID::prefitL15WithSatClock;
        lc15PrefitOfGALWithSatClock.body[TypeID::LC15]          =   +1.0;
        lc15PrefitOfGALWithSatClock.body[TypeID::rho]           =   -1.0;
        lc15PrefitOfGALWithSatClock.body[TypeID::cdtSat]        =   +1.0;
        lc15PrefitOfGALWithSatClock.body[TypeID::relativity]    =   -1.0;
        lc15PrefitOfGALWithSatClock.body[TypeID::gravDelay]     =   -1.0;
        lc15PrefitOfGALWithSatClock.body[TypeID::tropoSlant]    =   -1.0;
        lc15PrefitOfGALWithSatClock.body[TypeID::windUp]        =   -LC_WAVELENGTH_GAL_L1L5/TWO_PI;

         // Definition to compute prefit residual of Galileo LC, L1 + L5, with sta clock
        lc15PrefitOfGALWithStaClock.header                      =   TypeID::prefitL15WithStaClock;
        lc15PrefitOfGALWithStaClock.body[TypeID::LC]            =   +1.0;
        lc15PrefitOfGALWithStaClock.body[TypeID::rho]           =   -1.0;
        lc15PrefitOfGALWithStaClock.body[TypeID::cdtSta]        =   -1.0;
        lc15PrefitOfGALWithStaClock.body[TypeID::relativity]    =   -1.0;
        lc15PrefitOfGALWithStaClock.body[TypeID::gravDelay]     =   -1.0;
        lc15PrefitOfGALWithStaClock.body[TypeID::tropoSlant]    =   -1.0;
        lc15PrefitOfGALWithStaClock.body[TypeID::windUp]        =   -LC_WAVELENGTH_GAL_L1L5/TWO_PI;

        // Definition to compute prefit residual of Galileo LC, L1 + L5, for PCE
        lc15PrefitOfGALForPCE.header                            =   TypeID::prefitL15ForPCE;
        lc15PrefitOfGALForPCE.body[TypeID::LC15]                =   +1.0;
        lc15PrefitOfGALForPCE.body[TypeID::rho]                 =   -1.0;
        lc15PrefitOfGALForPCE.body[TypeID::relativity]          =   -1.0;
        lc15PrefitOfGALForPCE.body[TypeID::gravDelay]           =   -1.0;
        lc15PrefitOfGALForPCE.body[TypeID::tropoSlant]          =   -1.0;
        lc15PrefitOfGALForPCE.body[TypeID::windUp]              =   -LC_WAVELENGTH_GAL_L1L5/TWO_PI;

        // Definition to compute prefit residual of Galileo LC, L1 + L5, for POD
        lc15PrefitOfGALForPOD.header                            =   TypeID::prefitL15ForPOD;
        lc15PrefitOfGALForPOD.body[TypeID::LC15]                =   +1.0;
        lc15PrefitOfGALForPOD.body[TypeID::relativity]          =   -1.0;
        lc15PrefitOfGALForPOD.body[TypeID::gravDelay]           =   -1.0;
        lc15PrefitOfGALForPOD.body[TypeID::tropoSlant]          =   -1.0;
        lc15PrefitOfGALForPOD.body[TypeID::windUp]              =   -LC_WAVELENGTH_GAL_L1L5/TWO_PI;

        // Definition to compute prefit residual of Galileo LC, L1 + L5, for PPP
        lc15PrefitOfGALForPPP.header                            =   TypeID::prefitLForPPP;
        lc15PrefitOfGALForPPP.body[TypeID::LC15]                =   +1.0;
        lc15PrefitOfGALForPPP.body[TypeID::rho]                 =   -1.0;
        lc15PrefitOfGALForPPP.body[TypeID::cdtSat]              =   +1.0;
        lc15PrefitOfGALForPPP.body[TypeID::relativity]          =   -1.0;
        lc15PrefitOfGALForPPP.body[TypeID::gravDelay]           =   -1.0;
        lc15PrefitOfGALForPPP.body[TypeID::tropoSlant]          =   -1.0;
        lc15PrefitOfGALForPPP.body[TypeID::windUp]              =   -LC_WAVELENGTH_GAL_L1L5/TWO_PI;




         // Definition to compute prefit residual of BDS LC, L2 + L7
        lc27PrefitOfBDS.header                                  =   TypeID::prefitL27;
        lc27PrefitOfBDS.body[TypeID::LC27]                      =   +1.0;
        lc27PrefitOfBDS.body[TypeID::rho]                       =   -1.0;
        lc27PrefitOfBDS.body[TypeID::cdtSat]                    =   +1.0;
        lc27PrefitOfBDS.body[TypeID::relativity]                =   -1.0;
        lc27PrefitOfBDS.body[TypeID::gravDelay]                 =   -1.0;
        lc27PrefitOfBDS.body[TypeID::tropoSlant]                =   -1.0;
        lc27PrefitOfBDS.body[TypeID::windUp]                    =   -LC_WAVELENGTH_BDS_L2L7/TWO_PI;

         // Definition to compute prefit residual of BDS LC, L2 + L7, without clock
        lc27PrefitOfBDSWithoutClock.header                      =   TypeID::prefitL27WithoutClock;
        lc27PrefitOfBDSWithoutClock.body[TypeID::LC27]          =   +1.0;
        lc27PrefitOfBDSWithoutClock.body[TypeID::rho]           =   -1.0;
        lc27PrefitOfBDSWithoutClock.body[TypeID::relativity]    =   -1.0;
        lc27PrefitOfBDSWithoutClock.body[TypeID::gravDelay]     =   -1.0;
        lc27PrefitOfBDSWithoutClock.body[TypeID::tropoSlant]    =   -1.0;
        lc27PrefitOfBDSWithoutClock.body[TypeID::windUp]        =   -LC_WAVELENGTH_BDS_L2L7/TWO_PI;

         // Definition to compute prefit residual of BDS LC, L2 + L7, with sat clock
        lc27PrefitOfBDSWithSatClock.header                      =   TypeID::prefitL27WithSatClock;
        lc27PrefitOfBDSWithSatClock.body[TypeID::LC27]          =   +1.0;
        lc27PrefitOfBDSWithSatClock.body[TypeID::rho]           =   -1.0;
        lc27PrefitOfBDSWithSatClock.body[TypeID::cdtSat]        =   +1.0;
        lc27PrefitOfBDSWithSatClock.body[TypeID::relativity]    =   -1.0;
        lc27PrefitOfBDSWithSatClock.body[TypeID::gravDelay]     =   -1.0;
        lc27PrefitOfBDSWithSatClock.body[TypeID::tropoSlant]    =   -1.0;
        lc27PrefitOfBDSWithSatClock.body[TypeID::windUp]        =   -LC_WAVELENGTH_BDS_L2L7/TWO_PI;

         // Definition to compute prefit residual of BDS LC, L2 + L7, with sta clock
        lc27PrefitOfBDSWithStaClock.header                      =   TypeID::prefitL27WithStaClock;
        lc27PrefitOfBDSWithStaClock.body[TypeID::LC27]          =   +1.0;
        lc27PrefitOfBDSWithStaClock.body[TypeID::rho]           =   -1.0;
        lc27PrefitOfBDSWithStaClock.body[TypeID::cdtSta]        =   -1.0;
        lc27PrefitOfBDSWithStaClock.body[TypeID::relativity]    =   -1.0;
        lc27PrefitOfBDSWithStaClock.body[TypeID::gravDelay]     =   -1.0;
        lc27PrefitOfBDSWithStaClock.body[TypeID::tropoSlant]    =   -1.0;
        lc27PrefitOfBDSWithStaClock.body[TypeID::windUp]        =   -LC_WAVELENGTH_BDS_L2L7/TWO_PI;

        // Definition to compute prefit residual of BDS LC, L2 + L7, for PCE
        lc27PrefitOfBDSForPCE.header                            =   TypeID::prefitL27ForPCE;
        lc27PrefitOfBDSForPCE.body[TypeID::LC27]                =   +1.0;
        lc27PrefitOfBDSForPCE.body[TypeID::relativity]          =   -1.0;
        lc27PrefitOfBDSForPCE.body[TypeID::gravDelay]           =   -1.0;
        lc27PrefitOfBDSForPCE.body[TypeID::tropoSlant]          =   -1.0;
        lc27PrefitOfBDSForPCE.body[TypeID::windUp]              =   -LC_WAVELENGTH_BDS_L2L7/TWO_PI;

        // Definition to compute prefit residual of BDS LC, L2 + L7, for POD
        lc27PrefitOfBDSForPOD.header                            =   TypeID::prefitL27ForPOD;
        lc27PrefitOfBDSForPOD.body[TypeID::LC27]                =   +1.0;
        lc27PrefitOfBDSForPOD.body[TypeID::relativity]          =   -1.0;
        lc27PrefitOfBDSForPOD.body[TypeID::gravDelay]           =   -1.0;
        lc27PrefitOfBDSForPOD.body[TypeID::tropoSlant]          =   -1.0;
        lc27PrefitOfBDSForPOD.body[TypeID::windUp]              =   -LC_WAVELENGTH_BDS_L2L7/TWO_PI;

        // Definition to compute prefit residual of BDS LC, L2 + L7, for PPP
        lc27PrefitOfBDSForPPP.header                            =   TypeID::prefitL27ForPPP;
        lc27PrefitOfBDSForPPP.body[TypeID::LC27]                =   +1.0;
        lc27PrefitOfBDSForPPP.body[TypeID::cdtSat]              =   +1.0;
        lc27PrefitOfBDSForPPP.body[TypeID::relativity]          =   -1.0;
        lc27PrefitOfBDSForPPP.body[TypeID::gravDelay]           =   -1.0;
        lc27PrefitOfBDSForPPP.body[TypeID::tropoSlant]          =   -1.0;
        lc27PrefitOfBDSForPPP.body[TypeID::windUp]              =   -LC_WAVELENGTH_BDS_L2L7/TWO_PI;


         // Definition to compute prefit residual of GPS MW, L1+L2
        mw12PrefitOfGPS.header                              =   TypeID::prefitMW12;
        mw12PrefitOfGPS.body[TypeID::MW12]                  =   +1.0;
        mw12PrefitOfGPS.body[TypeID::instMW12]              =   +1.0;

         // Definition to compute prefit residual of GPS MW, L1+L5
        mw15PrefitOfGPS.header                              =   TypeID::prefitMW15;
        mw15PrefitOfGPS.body[TypeID::MW15]                  =   +1.0;
        mw15PrefitOfGPS.body[TypeID::instMW15]              =   +1.0;

         // Definition to compute prefit residual of Galileo MW, L1+L5
        mw15PrefitOfGAL.header                              =   TypeID::prefitMW15;
        mw15PrefitOfGAL.body[TypeID::MW15]                  =   +1.0;
        mw15PrefitOfGAL.body[TypeID::instMW15]              =   +1.0;

         // Definition to compute prefit residual of BDS MW, L2+L7
        mw27PrefitOfBDS.header                              =   TypeID::prefitMW27;
        mw27PrefitOfBDS.body[TypeID::MW27]                  =   +1.0;
        mw27PrefitOfBDS.body[TypeID::instMW27]              =   +1.0;


    }  // End of constructor 'LinearCombinations::LinearCombinations()'


      // Return the frequency of the combination in cycles: i * L1 + j * L2
    double LinearCombinations::freqOfLC(int i, int j, double f1 , double f2 )
    {
        return ( double(i)*f1+double(j)*f2 );
    }

    /// Return the wavelength of the combination in cycles: i * L1 + j * L2
    double LinearCombinations::wavelengthOfLC(int i,int j,double f1,double f2)
    {
        return C_MPS / freqOfLC(i,j,f1,f2);
    }

    /// Return the f1 factor of the combination in cycles: i * L1 + j * L2
    double LinearCombinations::firstFactorOfLC(int i,int j,double f1,double f2)
    {
        return double(i)*f1/freqOfLC(i,j,f1,f2);
    }

    /// Return the f2 factor of the combination in cycles: i * L1 + j * L2
    double LinearCombinations::secondFactorOfLC(int i,int j,double f1,double f2 )
    {
        return double(j)*f2/freqOfLC(i,j,f1,f2);
    }

} // End of namespace gpstk
