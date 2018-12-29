#pragma ident "$Id$"

/**
 * @file ComputePartials.cpp
 */

#include "ComputePartials.hpp"
#include "constants.hpp"


using namespace std;


namespace gpstk
{

    // Return a string identifying this object.
    std::string ComputePartials::getClassName() const
    { return "ComputePartials"; }


      /* Return a satTypeValueMap object, adding the new data generated when
       * calling a modeling object.
       *
       * @param time      Epoch.
       * @param gData     Data object holding the data.
       */
    satTypeValueMap& ComputePartials::Process( const CommonTime& time,
                                               satTypeValueMap& gData )
        throw(ProcessingException)
    {

        try
        {
            CommonTime utc( pRefSys->GPS2UTC(time) );

            SatIDSet satRejectedSet;

            double rho(0.0), rel(0.0);

            // satellite state
            Vector<double> state;

            Vector<double> rsat(3,0.0), vsat(3,0.0);

            // satellite pos/vel in ICRS
            posSatECI.resize(3,0.0);
            velSatECI.resize(3,0.0);

            // satellite pos/vel in ITRS
            posSatECF.resize(3,0.0);
            velSatECF.resize(3,0.0);

            double dclock( cdtSource/C_MPS );

            Matrix<double> c2tRaw, c2tDot;

            // Loop through all the satellites
            for(satTypeValueMap::iterator it = gData.begin();
                it != gData.end();
                ++it)
            {
                SatID sat( it->first );

//                if(satStateMap.find(sat) == satStateMap.end())
//                {
//                    satRejectedSet.insert( sat );
//                    continue;
//                }

//                // satellite state
//                state = satStateMap[sat];
//
//                // rsat_c, vsat_c
//                for(int i=0; i<3; ++i) rsat(i) = state(i+0);
//                for(int i=0; i<3; ++i) vsat(i) = state(i+3);
//
//                // satellite pos at (t-dT) in ICRS
//                posSatECI = rsat;
//
//                // satellite vel at (t-dT) in ICRS
//                velSatECI = vsat;
//
//                // norm( rsat(t) )
//                double r( norm(rsat) );
//                double r2( r * r );
//                double r3( r2 * r );
//                double k( GM_EARTH/r3 );
//
//                // initial value of rho
//                double rho( norm(posSatECI - posSourceECI) );
//                double dlight( rho/C_MPS );
//
//                double drel( 0.0 );
//
//                // rho/c + dtr + drel
//                double dT(dlight + dclock - drel);
//                double dT2(dT * dT);
//
//                // iterate to get rsat(t-dT)/vsat(t-dT)
//                for(int i=0; i<3; ++i)
//                {
//                    velSatECI = vsat + k*rsat*dT;
//                    posSatECI = rsat - vsat*dT - 0.5*k*rsat*dT2;
//
//                    rho = norm(posSatECI - posSourceECI);
//                    dlight = rho/C_MPS;
//
//                    drel = -2.0*( posSatECI(0)*velSatECI(0) +
//                                  posSatECI(1)*velSatECI(1) +
//                                  posSatECI(2)*velSatECI(2) )/C_MPS/C_MPS;
//
//                    dT = dlight + dclock - drel;
//                    dT2 = dT * dT;
//                }
//
//                // R(t-dT), dR(t-dT)
//                c2tRaw = pRefSys->C2TMatrix(utc-dT);
//                c2tDot = pRefSys->dC2TMatrix(utc-dT);
//
//                // satellite pos at (t-dT) in ITRS
//                posSatECF = c2tRaw * posSatECI;
//
//                // satellite vel at (t-dT) in ITRS
//                velSatECF = c2tRaw * velSatECI + c2tDot * posSatECI;
//
//                Position posSat( posSatECF(0), posSatECF(1), posSatECF(2) );
//
//                if(satPosECIMap.find(sat) == satPosECIMap.end())
//                    satPosECIMap[sat] = posSatECI;
//
//                if(satPosECFMap.find(sat) == satPosECFMap.end())
//                    satPosECFMap[sat] = posSatECF;
//
//                Matrix<double> I = ident<double>(3);
//
//                // rsat_c(t) * rsat_c(t)'
//                Matrix<double> rrt( outer(rsat,rsat) );
//
//                // drsat_c(t-dT) / drsat_c(t)
//                Matrix<double> phi( I - 0.5*k*(I-3.0*rrt/r2)*dT2 );
//
//                // rho
//                rho = norm(posSatECI - posSourceECI);
//
//                // drho / drsat_c(t-dT)
//                Vector<double> par( (posSatECI-posSourceECI)/rho );
//
//                // drho / drsat_c(t)
//                Vector<double> rho_rsat( phi * par );
//
//                // drho / dvsat_c(t)
//                Vector<double> rho_vsat( -par * dT );
//
//                // satellite position in ITRS
//                it->second[TypeID::satXECF] = posSatECF[0];
//                it->second[TypeID::satYECF] = posSatECF[1];
//                it->second[TypeID::satZECF] = posSatECF[2];
//
//                // satellite position in ICRS
//                it->second[TypeID::satXECI] = posSatECI[0];
//                it->second[TypeID::satYECI] = posSatECI[1];
//                it->second[TypeID::satZECI] = posSatECI[2];
//
//                // satellite velocity in ITRS
//                it->second[TypeID::satVXECF] = velSatECF[0];
//                it->second[TypeID::satVYECF] = velSatECF[1];
//                it->second[TypeID::satVZECF] = velSatECF[2];
//
//                // satellite velocity in ICRS
//                it->second[TypeID::satVXECI] = velSatECI[0];
//                it->second[TypeID::satVYECI] = velSatECI[1];
//                it->second[TypeID::satVZECI] = velSatECI[2];
//
//                // coefficients of satellite position in ICRS
//                it->second[TypeID::dSatX] = rho_rsat[0];
//                it->second[TypeID::dSatY] = rho_rsat[1];
//                it->second[TypeID::dSatZ] = rho_rsat[2];
//
//                // coefficients of satellite velocity in ICRS
//                it->second[TypeID::dSatVX] = rho_vsat[0];
//                it->second[TypeID::dSatVY] = rho_vsat[1];
//                it->second[TypeID::dSatVZ] = rho_vsat[2];


            } // End of loop for(stv = gData.begin()...

            // Remove satellites with missing data
            gData.removeSatID(satRejectedSet);

            return gData;

        }   // End of try...
        catch(Exception& u)
        {
            // Throw an exception if something unexpected happens
            ProcessingException e( getClassName() + ":" + u.what() );
            GPSTK_THROW(e);
        }

    }  // End of method 'ComputePartials::Process()'


     /** Return a gnssDataMap object, adding the new data generated when
      *  calling a modeling object.
      *
      * @param gData    Data object holding the data.
      */
    gnssDataMap& ComputePartials::Process(gnssDataMap& gData)
        throw(ProcessingException)
    {
        SourceIDSet sourceRejectedSet;

        sourcePosECFMap.clear();
        sourcePosECIMap.clear();

        if(pRefSys == NULL) return gData;
        if(pMSCStore == NULL) return gData;
        if(sourceClockMap.empty()) return gData;

        Matrix<double> t2cRaw(3,3,0.0);

        for(gnssDataMap::iterator gdmIt = gData.begin();
            gdmIt != gData.end();
            ++gdmIt)
        {
            CommonTime gps( gdmIt->first );
            CommonTime utc( pRefSys->GPS2UTC(gps) );

            for(sourceDataMap::iterator sdmIt = gdmIt->second.begin();
                sdmIt != gdmIt->second.end();
                ++sdmIt)
            {
                posSourceECI.resize(3,0.0);
                posSourceECF.resize(3,0.0);

                SourceID source( sdmIt->first );
                string station( source.sourceName );

                MSCData mscData;

                gps.setTimeSystem( TimeSystem::Unknown );
                try
                {
                    mscData = pMSCStore->findMSC(station,gps);
                }
                catch(...)
                {
                    sourceRejectedSet.insert(source);
                    continue;
                }
                gps.setTimeSystem( TimeSystem::GPS );

                nominalPos = mscData.coordinates;

                if(sourceClockMap.find(source) == sourceClockMap.end())
                {
                    sourceRejectedSet.insert(source);
                    continue;
                }

                posSourceECF[0] = nominalPos.X();
                posSourceECF[1] = nominalPos.Y();
                posSourceECF[2] = nominalPos.Z();

                cdtSource = sourceClockMap[source];

                t2cRaw = pRefSys->T2CMatrix(utc-cdtSource/C_MPS);

                posSourceECI = t2cRaw * posSourceECF;

                Process( gdmIt->first, sdmIt->second );

                sourcePosECFMap[source] = posSourceECF;
                sourcePosECIMap[source] = posSourceECI;
            }
        }

        gData.removeSourceID( sourceRejectedSet );

        return gData;

    }  // End of method 'ComputePartials::Process()'


}  // End of namespace gpstk
