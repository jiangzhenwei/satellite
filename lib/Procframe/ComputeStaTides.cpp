#pragma ident "$Id$"

/**
 * @file ComputeStaTides.cpp
 * This class computes station tides.
 */


#include "ComputeStaTides.hpp"

using namespace std;

namespace gpstk
{

    // Return a string identifying this object.
    std::string ComputeStaTides::getClassName() const
    { return "ComputeStaTides"; }



    Triple ComputeStaTides::getTides( const CommonTime& time,
                                      const string& station )
        throw(ProcessingException)
    {

        try
        {
            Triple solidCorr, oceanCorr, poleCorr;

            solidCorr = solid.getSolidTide(time, nominalPos);
            oceanCorr = ocean.getOceanLoading(time, station);
            poleCorr = pole.getPoleTide(time, nominalPos);

            return solidCorr + oceanCorr + poleCorr;
        }
        catch(Exception& u)
        {
            // Throw an exception if something unexpected happens
            ProcessingException e( getClassName() + ":" + u.what() );
            GPSTK_THROW(e);
        }

    }  // End of method 'ComputeStaTides::Process()'


}  // End of namespace gpstk
