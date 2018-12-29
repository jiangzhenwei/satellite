#pragma ident "$Id$"

/**
 * @file LinearCombinations.hpp
 * This class defines handy linear combinations of GDS data.
 */

#ifndef LINEARCOMBINATIONS_HPP
#define LINEARCOMBINATIONS_HPP

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


#include "DataStructures.hpp"
#include "constants.hpp"


namespace gpstk
{

      /** @addtogroup DataStructures */
      //@{


      /** This class defines handy linear combinations of GDS data.
       *
       * This class is meant to be used with the GNSS data structures (GDS)
       * objects found in "DataStructures" class, and it is intended to be
       * coupled with class ComputeLinear.hpp.
       *
       * A typical way to use this class follows:
       *
       * @code
       *
       *      // Define a LinearCombinations object
       *   LinearCombinations comb;
       *
       *      // Object to compute linear combinations of data
       *      // Linear combinations will be computed in a FIFO basis
       *   ComputeLinear linear;
       *
       *      // Add a linear combination to compute PC combination using C1
       *   linear.addLinear(comb.pcCombWithC1);
       *
       *      // Add a linear combination to compute prefit residual using PC
       *   linear.addLinear(comb.pcPrefit);
       *
       *
       *      // Load observation data
       *   RinexObsStream rin("ebre0300.02o");
       *
       *      // Loads precise ephemeris object with file data
       *   SP3EphemerisStore SP3EphList;
       *   SP3EphList.loadFile("igs11513.sp3");
       *
       *      // Sets nominal position of receiver
       *   Position nominalPos(4833520.3800, 41536.8300, 4147461.2800);
       *
       *      // Declare a MOPSTropModel object, setting the defaults
       *   MOPSTropModel mopsTM( nominalPos.getAltitude(),
       *                         nominalPos.getGeodeticLatitude(), 30);
       *
       *      // Object to compute the tropospheric data
       *   ComputeTropModel computeTropo(mopsTM);
       *
       *      // Declare a basic modeler
       *   BasicModel basic(nominalPos, SP3EphList);
       *
       *   gnssRinex gRin;
       *
       *   while(rin >> gRin)
       *   {
       *
       *      gRin >> basic >> computeTropo >> linear;
       *
       *         // Dump results
       *      gRin.body.dump(cout,1);
       *
       *   }
       *
       * @endcode
       *
       * @sa ComputeLinear.hpp
       */
    class LinearCombinations
    {
    public:

        /// Default constructor
        LinearCombinations();



        /// Definition to compute MIGC1 (GPS L1 with minus ionospehric delays)
        gnssLinearCombination q1CombOfGPSL1L2;

        /// Definition to compute MIGC2 (GPS L2 with minus ionospehric delays)
        gnssLinearCombination q2CombOfGPSL1L2;


        /// Definition to compute MIGC1 (GPS L1 with minus ionospehric delays)
        gnssLinearCombination q1CombOfGPSL1L5;

        /// Definition to compute MIGC5 (GPS L5 with minus ionospehric delays)
        gnssLinearCombination q5CombOfGPSL1L5;


        /// Definition to compute MIEC1 (Galileo L1 with minus ionospehric delays)
        gnssLinearCombination q1CombOfGALL1L5;

        /// Definition to compute MIEC5 (Galileo L5 with minus ionospehric delays)
        gnssLinearCombination q5CombOfGALL1L5;


        /// Definition to compute MICC2 (BDS L2 with minus ionospehric delays)
        gnssLinearCombination q2CombOfBDSL2L7;

        /// Definition to compute MICC7 (BDS L7 with minus ionospehric delays)
        gnssLinearCombination q7CombOfBDSL2L7;



        /// Definition to compute PC combination of GPS, C1 and C2
        gnssLinearCombination pc12CombOfGPS;

        /// Definition to compute LC combination of GPS, L1 and L2
        gnssLinearCombination lc12CombOfGPS;


        /// Definition to compute PC combination of GPS, C1 and C5
        gnssLinearCombination pc15CombOfGPS;

        /// Definition to compute LC combination of GPS, L1 and L5
        gnssLinearCombination lc15CombOfGPS;


        /// Definition to compute PC combination of Galileo, C1 and C5
        gnssLinearCombination pc15CombOfGAL;

        /// Definition to compute LC combination of Galileo, L1 and L5
        gnssLinearCombination lc15CombOfGAL;


        /// Definition to compute PC combination of BDS, C2 and C7
        gnssLinearCombination pc27CombOfBDS;

        /// Definition to compute LC combination of BDS, L2 and L7
        gnssLinearCombination lc27CombOfBDS;



        /// Definition to compute PI combination of GPS, C1 and C2
        gnssLinearCombination pi12CombOfGPS;

        /// Definition to compute LI combination of GPS, L1 and L2
        gnssLinearCombination li12CombOfGPS;


        /// Definition to compute PI combination of GPS, C1 and C5
        gnssLinearCombination pi15CombOfGPS;

        /// Definition to compute LI combination of GPS, L1 and L5
        gnssLinearCombination li15CombOfGPS;


        /// Definition to compute PI combination of Galileo, C1 and C5
        gnssLinearCombination pi15CombOfGAL;

        /// Definition to compute LI combination of Galileo, L1 and L5
        gnssLinearCombination li15CombOfGAL;


        /// Definition to compute PI combination of BDS, C2 and C7
        gnssLinearCombination pi27CombOfBDS;

        /// Definition to compute LI combination of BDS, L2 and L7
        gnssLinearCombination li27CombOfBDS;



        /// Definition to compute PW combination of GPS, C1 and C2
        gnssLinearCombination pw12CombOfGPS;

        /// Definition to compute LW combination of GPS, L1 and L2
        gnssLinearCombination lw12CombOfGPS;

        /// Definition to compute PW combination of GPS, C1 and C5
        gnssLinearCombination pw15CombOfGPS;

        /// Definition to compute LW combination of GPS, L1 and L5
        gnssLinearCombination lw15CombOfGPS;


        /// Definition to compute PW combination of Galileo, C1 and C5
        gnssLinearCombination pw15CombOfGAL;

        /// Definition to compute LW combination of Galileo, L1 and L5
        gnssLinearCombination lw15CombOfGAL;


        /// Definition to compute PW combination of BDS, C2 and C7
        gnssLinearCombination pw27CombOfBDS;

        /// Definition to compute LW combination of BDS, L2 and L7
        gnssLinearCombination lw27CombOfBDS;



        /// Definition to compute PN combination of GPS, C1 and C2
        gnssLinearCombination pn12CombOfGPS;

        /// Definition to compute LN combination of GPS, L1 and L2
        gnssLinearCombination ln12CombOfGPS;


        /// Definition to compute PN combination of GPS, C1 and C5
        gnssLinearCombination pn15CombOfGPS;

        /// Definition to compute LN combination of GPS, L1 and L5
        gnssLinearCombination ln15CombOfGPS;


        /// Definition to compute PN combination of Galileo, C1 and C5
        gnssLinearCombination pn15CombOfGAL;

        /// Definition to compute LN combination of Galileo, L1 and L5
        gnssLinearCombination ln15CombOfGAL;


        /// Definition to compute PN combination of BDS, C2 and C7
        gnssLinearCombination pn27CombOfBDS;

        /// Definition to compute LN combination of BDS, L2 and L7
        gnssLinearCombination ln27CombOfBDS;



        /// Definition to compute the MW combination of GPS, C1+C2 and L1+L2
        gnssLinearCombination mw12CombOfGPS;

        /// Definition to compute the MW combination of GPS, C1+C5 and L1+L5
        gnssLinearCombination mw15CombOfGPS;

        /// Definition to compute the MW combination of Galileo, C1+C5 and L1+L5
        gnssLinearCombination mw15CombOfGAL;

        /// Definition to compute the MW combination of BDS, C2+C7 and L2+L7
        gnssLinearCombination mw27CombOfBDS;



        /// Definition to compute prefit residual of GPS
        gnssLinearCombination c1PrefitOfGPS;
        gnssLinearCombination c2PrefitOfGPS;
        gnssLinearCombination c5PrefitOfGPS;


        /// Definition to compute prefit residual of Galileo
        gnssLinearCombination c1PrefitOfGAL;
        gnssLinearCombination c5PrefitOfGAL;
        gnssLinearCombination c6PrefitOfGAL;
        gnssLinearCombination c7PrefitOfGAL;
        gnssLinearCombination c8PrefitOfGAL;


        /// Definition to compute prefit residual of BDS
        gnssLinearCombination c2PrefitOfBDS;
        gnssLinearCombination c6PrefitOfBDS;
        gnssLinearCombination c7PrefitOfBDS;



        /// Definition to compute prefit residual of GPS
        gnssLinearCombination l1PrefitOfGPS;
        gnssLinearCombination l2PrefitOfGPS;
        gnssLinearCombination l5PrefitOfGPS;


        /// Definition to compute prefit residual of Galileo
        gnssLinearCombination l1PrefitOfGAL;
        gnssLinearCombination l5PrefitOfGAL;
        gnssLinearCombination l6PrefitOfGAL;
        gnssLinearCombination l7PrefitOfGAL;
        gnssLinearCombination l8PrefitOfGAL;


        /// Definition to compute prefit residual of BDS
        gnssLinearCombination l2PrefitOfBDS;
        gnssLinearCombination l6PrefitOfBDS;
        gnssLinearCombination l7PrefitOfBDS;


        /// Definition to compute prefit residual of GPS
        gnssLinearCombination pc12PrefitOfGPS;
        gnssLinearCombination pc12PrefitOfGPSWithoutClock;
        gnssLinearCombination pc12PrefitOfGPSWithSatClock;
        gnssLinearCombination pc12PrefitOfGPSWithStaClock;
        gnssLinearCombination pc12PrefitOfGPSForPCE;
        gnssLinearCombination pc12PrefitOfGPSForPOD;
        gnssLinearCombination pc12PrefitOfGPSForPPP;

        gnssLinearCombination lc12PrefitOfGPS;
        gnssLinearCombination lc12PrefitOfGPSWithoutClock;
        gnssLinearCombination lc12PrefitOfGPSWithSatClock;
        gnssLinearCombination lc12PrefitOfGPSWithStaClock;
        gnssLinearCombination lc12PrefitOfGPSForPCE;
        gnssLinearCombination lc12PrefitOfGPSForPOD;
        gnssLinearCombination lc12PrefitOfGPSForPPP;

        gnssLinearCombination pc15PrefitOfGPS;
        gnssLinearCombination pc15PrefitOfGPSWithoutClock;
        gnssLinearCombination pc15PrefitOfGPSWithSatClock;
        gnssLinearCombination pc15PrefitOfGPSWithStaClock;
        gnssLinearCombination pc15PrefitOfGPSForPCE;
        gnssLinearCombination pc15PrefitOfGPSForPOD;
        gnssLinearCombination pc15PrefitOfGPSForPPP;

        gnssLinearCombination lc15PrefitOfGPS;
        gnssLinearCombination lc15PrefitOfGPSWithoutClock;
        gnssLinearCombination lc15PrefitOfGPSWithSatClock;
        gnssLinearCombination lc15PrefitOfGPSWithStaClock;
        gnssLinearCombination lc15PrefitOfGPSForPCE;
        gnssLinearCombination lc15PrefitOfGPSForPOD;
        gnssLinearCombination lc15PrefitOfGPSForPPP;


        /// Definition to compute prefit residual of Galileo
        gnssLinearCombination pc15PrefitOfGAL;
        gnssLinearCombination pc15PrefitOfGALWithoutClock;
        gnssLinearCombination pc15PrefitOfGALWithSatClock;
        gnssLinearCombination pc15PrefitOfGALWithStaClock;
        gnssLinearCombination pc15PrefitOfGALForPCE;
        gnssLinearCombination pc15PrefitOfGALForPOD;
        gnssLinearCombination pc15PrefitOfGALForPPP;

        gnssLinearCombination lc15PrefitOfGAL;
        gnssLinearCombination lc15PrefitOfGALWithoutClock;
        gnssLinearCombination lc15PrefitOfGALWithSatClock;
        gnssLinearCombination lc15PrefitOfGALWithStaClock;
        gnssLinearCombination lc15PrefitOfGALForPCE;
        gnssLinearCombination lc15PrefitOfGALForPOD;
        gnssLinearCombination lc15PrefitOfGALForPPP;


        /// Definition to compute prefit residual of BDS
        gnssLinearCombination pc27PrefitOfBDS;
        gnssLinearCombination pc27PrefitOfBDSWithoutClock;
        gnssLinearCombination pc27PrefitOfBDSWithSatClock;
        gnssLinearCombination pc27PrefitOfBDSWithStaClock;
        gnssLinearCombination pc27PrefitOfBDSForPCE;
        gnssLinearCombination pc27PrefitOfBDSForPOD;
        gnssLinearCombination pc27PrefitOfBDSForPPP;

        gnssLinearCombination lc27PrefitOfBDS;
        gnssLinearCombination lc27PrefitOfBDSWithoutClock;
        gnssLinearCombination lc27PrefitOfBDSWithSatClock;
        gnssLinearCombination lc27PrefitOfBDSWithStaClock;
        gnssLinearCombination lc27PrefitOfBDSForPCE;
        gnssLinearCombination lc27PrefitOfBDSForPOD;
        gnssLinearCombination lc27PrefitOfBDSForPPP;


        /// Definition to compute prefit residual of GPS MW, L1+L2
        gnssLinearCombination mw12PrefitOfGPS;

        /// Definition to compute prefit residual of GPS MW, L1+L5
        gnssLinearCombination mw15PrefitOfGPS;


        /// Definition to compute prefit residual of Galileo MW, L1+L5
        gnssLinearCombination mw15PrefitOfGAL;

        /// Definition to compute prefit residual of BDS MW, L2+L7
        gnssLinearCombination mw27PrefitOfBDS;


    public:

        /// Return the frequency of the combination in cycles: i * L1 + j * L2
        static double freqOfLC(int i, int j, double f1 = L1_FREQ_GPS, double f2 = L2_FREQ_GPS);

        /// Return the wavelength of the combination in cycles: i * L1 + j * L2
        static double wavelengthOfLC(int i, int j, double f1 = L1_FREQ_GPS, double f2 = L2_FREQ_GPS);

        /// Return the f1 factor of the combination in cycles: i * L1 + j * L2
        static double firstFactorOfLC(int i, int j, double f1 = L1_FREQ_GPS, double f2 = L2_FREQ_GPS);

        /// Return the f2 factor of the combination in cycles: i * L1 + j * L2
        static double secondFactorOfLC(int i, int j, double f1 = L1_FREQ_GPS, double f2 = L2_FREQ_GPS);

    }; // End of class 'LinearCombinations'

    //@}

}

#endif // LINEARCOMBINATIONS_HPP
