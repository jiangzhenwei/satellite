#pragma ident "$Id$"

/**
 * @file ComputeDOP.cpp
 * This class computes the usual DOP values: GDOP, PDOP, TDOP, HDOP and VDOP.
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
//  Dagoberto Salazar - gAGE ( http://www.gage.es ). 2008, 2011
//
//============================================================================


#include "ComputeDOP.hpp"


namespace gpstk
{

      // Return a string identifying this object.
   std::string ComputeDOP::getClassName() const
   { return "ComputeDOP"; }



      /* Return a satTypeValueMap object, adding the new data generated when
       * calling this object.
       *
       * @param time      Epoch corresponding to the data.
       * @param gData     Data object holding the data.
       */
   satTypeValueMap& ComputeDOP::Process( const CommonTime& time,
                                         satTypeValueMap& gData)
      throw(ProcessingException)
   {

      try
      {

         bool valid1(false), valid2(false);

            // First, let's define a set with XYZt unknowns
         TypeIDSet tempSet1;
         tempSet1.insert(TypeID::dStaX);
         tempSet1.insert(TypeID::dStaY);
         tempSet1.insert(TypeID::dStaZ);
         tempSet1.insert(TypeID::cdtSta);

            // Second, let's define a set with NEUt unknowns
         TypeIDSet tempSet2;
         tempSet2.insert(TypeID::dStaLat);
         tempSet2.insert(TypeID::dStaLon);
         tempSet2.insert(TypeID::dStaH);
         tempSet2.insert(TypeID::cdtSta);

            // Then, generate the corresponding geometry/design matrices
         Matrix<double> dMatrix1(gData.getMatrixOfTypes(tempSet1));
         Matrix<double> dMatrix2(gData.getMatrixOfTypes(tempSet2));

            // Afterwards, compute the appropriate extra matrices
         Matrix<double> AT1(transpose(dMatrix1));
         Matrix<double> covM1(AT1 * dMatrix1);

         Matrix<double> AT2(transpose(dMatrix2));
         Matrix<double> covM2(AT2 * dMatrix2);

            // Let's try to invert AT*A matrices
         try
         {

            covM1 = inverseChol( covM1 );
            valid1 = true;

         }
         catch(...)
         {

            valid1 = false;
         }

         try
         {

            covM2 = inverseChol( covM2 );
            valid2 = true;

         }
         catch(...)
         {
            valid2 = false;
         }

         if( valid1 )
         {

            gdop = std::sqrt(covM1(0,0)+covM1(1,1)+covM1(2,2)+covM1(3,3));
            pdop = std::sqrt(covM1(0,0)+covM1(1,1)+covM1(2,2));
            tdop = std::sqrt(covM1(3,3));

         }
         else
         {
            gdop = -1.0;
            pdop = -1.0;
            tdop = -1.0;
         }

         if( valid2 )
         {
            hdop = std::sqrt(covM2(0,0)+covM2(1,1));
            vdop = std::sqrt(covM2(2,2));
         }
         else
         {
            hdop = -1.0;
            vdop = -1.0;
         }

         return gData;

      }
      catch(Exception& u)
      {
            // Throw an exception if something unexpected happens
         ProcessingException e( getClassName() + ":"
                                + u.what() );

         GPSTK_THROW(e);

      }

   }  // End of method 'ComputeDOP::Process()'


     /** Return a gnssDataMap object, adding the new data generated
      *  when calling this object.
      *
      * @param gData    Data object holding the data.
      */
   gnssDataMap& ComputeDOP::Process(gnssDataMap& gData)
      throw(ProcessingException)
   {

      for( gnssDataMap::iterator gdmIt = gData.begin();
           gdmIt != gData.end();
           ++gdmIt )
      {
         for( sourceDataMap::iterator sdmIt = gdmIt->second.begin();
              sdmIt != gdmIt->second.end(); sdmIt++ )
         {
            Process( gdmIt->first, sdmIt->second );
         }
      }

      return gData;

   }  // End of method 'ComputeDOP::Process()'

}  // End of namespace gpstk
