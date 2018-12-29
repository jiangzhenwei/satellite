#pragma ident "$Id$"

/**
 * @file ComputeStaTides.hpp
 * This class computes station tides.
 */

#ifndef COMPUTE_STA_TIDES_HPP
#define COMPUTE_STA_TIDES_HPP


#include <string>

#include "ProcessingClass.hpp"
#include "SolidTides.hpp"
#include "OceanLoading.hpp"
#include "PoleTides.hpp"


namespace gpstk
{

    class ComputeStaTides
    {
    public:

         /// Default constructor
        ComputeStaTides() {};


        Triple getTides( const CommonTime& time,
                         const std::string& station )
            throw(ProcessingException);


        virtual ComputeStaTides& setNominalPosition(const Position& pos)
        { nominalPos = pos; return (*this); };


        ComputeStaTides& setBLQDataReader( BLQDataReader& blqStore )
        { ocean.setBLQDataReader(blqStore); return (*this); };


        ComputeStaTides& setReferenceSystem( ReferenceSystem& refSys )
        {
            solid.setReferenceSystem(refSys);
            ocean.setReferenceSystem(refSys);
            pole.setReferenceSystem(refSys);
            return (*this);
        };


        ComputeStaTides& setSolarSystem( SolarSystem& solSys )
        {
            solid.setSolarSystem(solSys);
            return (*this);
        };


         /// Return a string identifying this object.
        virtual std::string getClassName(void) const;


         /// Destructor
        virtual ~ComputeStaTides() {};


    private:

        Position nominalPos;

        SolidTides solid;

        OceanLoading ocean;

        PoleTides pole;

    }; // End of class 'ComputeStaTides'

    //@}

}  // End of namespace gpstk

#endif  // COMPUTE_STA_TIDES_HPP
