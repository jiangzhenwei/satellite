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
//  Kaifa Kuang - Wuhan University . 2016
//
//============================================================================

/**
 * @file ECOM1Model.hpp
 * Class to do SRP calculation with ECOM1 model.
 */

#ifndef ECOM1_MODEL_HPP
#define ECOM1_MODEL_HPP

#include "SRPModel.hpp"


namespace gpstk
{
    /** @addtogroup GeoDynamics */
    //@{

    /** Class to do SRP calculation with ECOM1 model.
     */
    class ECOM1Model : public SRPModel
    {
    public:

        /// Default constructor
        ECOM1Model() {};

        /// Default destructor
        virtual ~ECOM1Model() {};


        /** Compute acceleration (and related partial derivatives) of SRP.
         * @param tt        TT
         * @param orbits    orbits
         */
        virtual void Compute( const CommonTime& tt,
                              const satVectorMap& orbits );

        /// Return the force model name
        inline virtual std::string modelName() const
        { return "ECOM1Model"; }

    }; // End of class 'ECOM1Model'

    // @}

}  // End of namespace 'gpstk'

#endif   // ECOM1_MODEL_HPP
