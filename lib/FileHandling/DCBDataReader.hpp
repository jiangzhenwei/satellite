#pragma ident "$Id$"

/**
 * @file DCBReader.hpp
 * Class to read DCB data from CODE.
 */

#ifndef DCBDATAREADER_HPP
#define DCBDATAREADER_HPP

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
//  Wei Yan - Chinese Academy of Sciences  2009, 2010
//
//============================================================================


#include <string>
#include <map>

#include "Exception.hpp"
#include "FFTextStream.hpp"
#include "StringUtils.hpp"
#include "SatID.hpp"


namespace gpstk
{

      /** @addtogroup formattedfile */
      //@{

      /** This is a class to read DCB (Differences of Code Biases) data
       *  from CODE.
       *
       * You can find DCB data at:
       *
       *    ftp.unibe.ch/aiub/BSWUSER50/ORB   - daily P1-P2
       *    ftp.unibe.ch/aiub/CODE            - monthly P1-P2 and P1-C1
       *
       *
       *  You should use different objects to load different DCB files. A typical
       *  way to use these classes follows:
       *
       * @code
       *      // Declare some DCBDataReader objects
       *   DCBDataReader dcbP1P2("P1P21002_ALL.DCB");
       *   DCBDataReader dcbP1C1("P1C11002.DCB");
       *
       *   double p1p2Sat1 = dcbP1P2.getDCB(1, SatID::systemGPS);
       *   double p1c1Sat1 = dcbP1C1.getDCB(1, SatID::systemGPS);
       *
       *   double p1p2ALGO = dcbP1P2.getDCB("ALGO");
       *
       * @endcode
       *
       * @sa DCBDataReader.hpp
       */
    class DCBDataReader : public FFTextStream
    {
    public:
         /// Default constructor
        DCBDataReader() {};

         /** Common constructor. It will always open file for read and will
          *  load DCB data in one pass.
          *
          * @param fn   DCB data file to read
          *
          */
        DCBDataReader(const char* fn)
            : FFTextStream(fn, std::ios::in)
        {
            FFTextStream(fn, std::ios::in);
            if( !FFTextStream::is_open() )
            {
                std::cerr << "The file " << fn << "doesn't exist!" << std::endl;
                exit(-1);
            }

            loadData();
        }


         /** Common constructor. It will always open file for read and will
          *  load DCB data in one pass.
          *
          * @param fn   DCB data file to read
          *
          */
        DCBDataReader(const std::string& fn)
        {
            FFTextStream(fn, std::ios::in);
            if( !FFTextStream::is_open() )
            {
                std::cerr << "The file " << fn << "doesn't exist!" << std::endl;
                exit(-1);
            }

            loadData();
        }

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Woverloaded-virtual"
        /// Method to open AND load DCB data file.
        virtual void open(const char* fn);


        /// Method to open AND load DCB data file. It doesn't
        /// clear data previously loaded.
        virtual void open(const std::string& fn);
#pragma clang diagnostic pop

        /// Get DCB data of a satellite.
        /// @param    sat   the satellite you desired
        /// @return         P1-P2 or P1-C1 depend what you have loaded
        double getDCB( const SatID& sat);

        /// Get DCB data of a satellite.
        /// @param    prn   the satellite id you desired
        /// @param    sys   the satellite system you desired
        /// @return         P1-P2 or P1-C1 depend what you have loaded
        double getDCB(const int& prn,
            const SatID::SatelliteSystem& sys = SatID::systemGPS);


        /// Get DCB data of a receiver.
        /// @param   sta   the receiver name you desired
        /// @param   sys   the satellite system you desired
        /// @return         P1-P2 or P1-C1 depend what you have loaded
        double getDCB( const std::string& sta,
                       const SatID::SatelliteSystem& sys = SatID::systemGPS);


        /// Destructor
        virtual ~DCBDataReader() {};


    private:

        // Map holding satellite DCB data
        typedef std::map< SatID, double > SatDCBData;

        // Map holding receiver DCB data
        typedef std::map< std::string, double > ReceiverDCBData;

        /// A structure used to store daily DCB data
        struct DailyDCBData
        {
            SatDCBData        satDCB;
            ReceiverDCBData   gpsDCB;
            ReceiverDCBData   glonassDCB;
        };

        /// Object holding all of the DCB data
        DailyDCBData allDCB;

        /// Method to store ocean tide harmonics data in this class' data map
        virtual void loadData()
            throw( FFStreamError, StringUtils::StringException );

    };  // End of class 'DCBDataReader'

    //@}

}  // End of namespace gpstk

#endif  // DCBDATAREADER_HPP
