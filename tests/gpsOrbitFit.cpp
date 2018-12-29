/* GNSS Orbit Fit. */

#include <iostream>

#include "ConfDataReader.hpp"

#include "DataStructures.hpp"

#include "ReferenceSystem.hpp"

#include "SolarSystem.hpp"

#include "SP3EphemerisStore.hpp"

#include "GNSSOrbit.hpp"

#include "RKF78Integrator.hpp"

#include "AdamsIntegrator.hpp"

#include "Epoch.hpp"

#include "StringUtils.hpp"

#include "Counter.hpp"


using namespace std;
using namespace gpstk;
using namespace gpstk::StringUtils;


int main(int argc, char* argv[])
{

    //---------- GLobal Data ----------//

    // Configuration File
    ConfDataReader confData;
    try
    {
        confData.open("orbit.conf");
    }
    catch(...)
    {
        cerr << "Configuration File Open Error." << endl;
        return 1;
    }


    // EOP File
    EOPDataStore2 eopDataStore;
    try
    {
        string eopFile = confData.getValue("eopFileName", "DEFAULT");
        eopDataStore.loadIERSFile(eopFile);
    }
    catch(...)
    {
        cerr << "EOP File Load Error." << endl;
        return 1;
    }

    eopDataStore.setInterpPoints(4);
    eopDataStore.setRegularization(true);


    // LeapSecond File
    LeapSecStore leapSecStore;
    try
    {
        string lsFile = confData.getValue("lsFileName", "DEFAULT");
        leapSecStore.loadFile(lsFile);
    }
    catch(...)
    {
        cerr << "LeapSecond File Load Error." << endl;
        return 1;
    }


    // Reference System
    ReferenceSystem refSys;
    refSys.setEOPDataStore(eopDataStore);
    refSys.setLeapSecStore(leapSecStore);


    // Solar System
    SolarSystem solSys;
    try
    {
        string ephFile = confData.getValue("jplFileName", "DEFAULT");
        solSys.initializeWithBinaryFile(ephFile);
    }
    catch(...)
    {
        cerr << "Solar System Initialize Error." << endl;
        return 1;
    }


    //---------- Force Model Configuration ----------//

    // Earth Gravitation
    EGM08Model egm(12,12);

    egm.setReferenceSystem(refSys);

    try
    {
        string egmFile = confData.getValue("egmFileName", "DEFAULT");
        egm.loadFile(egmFile);
    }
    catch(...)
    {
        cerr << "EGM File Load Error." << endl;
        return 1;
    }


    // Earth Solid Tide
    EarthSolidTide solidTide;
    solidTide.setReferenceSystem(refSys);
    solidTide.setSolarSystem(solSys);

    egm.setEarthSolidTide(solidTide);


    // Earth Ocean Tide
    EarthOceanTide oceanTide(8,8);
    oceanTide.setReferenceSystem(refSys);

    try
    {
        string eotFile = confData.getValue("otFileName", "DEFAULT");
        oceanTide.loadFile(eotFile);
    }
    catch(...)
    {
        cerr << "EOT File Load Error." << endl;
        return 1;
    }

    egm.setEarthOceanTide(oceanTide);


    // Earth Pole Tide
    EarthPoleTide poleTide;
    poleTide.setReferenceSystem(refSys);

    egm.setEarthPoleTide(poleTide);


    // Third Body
    ThirdBody thd;
    thd.setSolarSystem(solSys);
    thd.setReferenceSystem(refSys);
    thd.enableAllPlanets();


    // Solar Pressure
    int numSRP(5);
    Vector<double> srpc0(numSRP,0.0);

    satVectorMap satSRPC;

//    int maxNumSat(MAX_PRN_GPS + MAX_PRN_GAL + MAX_PRN_BDS);
    int maxNumSat(MAX_PRN_GPS);

    SatID sat;
    for(int i=1; i<=maxNumSat; ++i)
    {
        if(i <= MAX_PRN_GPS)
        {
            sat.id = i;
            sat.system = SatID::systemGPS;
        }
        else if(i <= MAX_PRN_GPS + MAX_PRN_GAL)
        {
            sat.id = i - MAX_PRN_GPS;
            sat.system = SatID::systemGalileo;
        }
        else
        {
            sat.id = i - MAX_PRN_GPS - MAX_PRN_GAL;
            sat.system = SatID::systemBDS;
        }

        satSRPC[sat] = srpc0;
    }

    ECOM1Model srp;
    srp.setReferenceSystem(refSys);
    srp.setSolarSystem(solSys);
    srp.setSRPCoeff(satSRPC);


    // Relativity
    Relativity rel;


    // GNSS Orbit
    GNSSOrbit gnss;
    gnss.setEGMModel(egm);
    gnss.setThirdBody(thd);
    gnss.setSRPModel(srp);
    gnss.setRelativity(rel);


    //---------- Integrator Configuration ----------//

    // RKF78 Integrator
    RKF78Integrator rkf78(60.0);

    rkf78.setEquationOfMotion(gnss);


    // Adams Integrator
    AdamsIntegrator adams(300.0);

    adams.setEquationOfMotion(gnss);


    // Arc Length
    double arcLen(24.0);

    // Arc Interval
    double arcInt(300.0);


    //---------- Initial Time ----------//
    string t0;
    try
    {
        t0 = confData.getValue("initialTime", "DEFAULT");
    }
    catch(...)
    {
        cerr << "Get Initial Time Error." << endl;
        return 1;
    }

    CommonTime gps0,utc0,tt0;

    CivilTime cv0;
    cv0.year    =    asInt( t0.substr( 0, 4) );
    cv0.month   =    asInt( t0.substr( 5, 2) );
    cv0.day     =    asInt( t0.substr( 8, 2) );
    cv0.hour    =    asInt( t0.substr(11, 2) );
    cv0.minute  =    asInt( t0.substr(14, 2) );
    cv0.second  = asDouble( t0.substr(17, 5) );

    cv0.setTimeSystem(TimeSystem::GPS);
    gps0 = cv0;
    utc0 = refSys.GPS2UTC(gps0);

    tt0 = refSys.GPS2TT(gps0);

    cout << "Initial Time: " << CivilTime(gps0) << endl;


    //---------- Reference Orbit ----------//

    // sp3 file
    string sp3FileListName;

    try
    {
        sp3FileListName = confData.getValue("sp3FileListName", "DEFAULT");
    }
    catch(...)
    {
        cerr << "sp3 file list name get error." << endl;
        exit(-1);
    }

    vector<string> sp3FileListVec;

    ifstream sp3FileListStream;

    sp3FileListStream.open(sp3FileListName.c_str());

    if(!sp3FileListStream)
    {
        cerr << "sp3 file list '" << sp3FileListName
             << "' open error." << endl;
        exit(-1);
    }

    string sp3File;
    while(sp3FileListStream >> sp3File)
    {
        sp3FileListVec.push_back(sp3File);
    }

    sp3FileListStream.close();

    if(sp3FileListVec.size() == 0)
    {
        cerr << "sp3 file list is empty." << endl;
        exit(-1);
    }

    SP3EphemerisStore sp3Store;
    sp3Store.rejectBadPositions(true);
    sp3Store.rejectBadClocks(true);
    sp3Store.setPosGapInterval(900+1);
    sp3Store.setPosMaxInterval(9*900+1);

    vector<string>::const_iterator sp3_it = sp3FileListVec.begin();
    while(sp3_it != sp3FileListVec.end())
    {
        sp3File = (*sp3_it);

        try
        {
            sp3Store.loadFile(sp3File);
        }
        catch(...)
        {
            cerr << "sp3 file '" << sp3File
                 << "' load error." << endl;

            ++sp3_it;
            continue;
        }

        ++sp3_it;
    }


    //---------- Orbit Fit ----------//

    Matrix<double> t2cRaw(3,3,0.0);
    Matrix<double> t2cDot(3,3,0.0);

    Vector<double> r0_itrs(3,0.0);
    Vector<double> v0_itrs(3,0.0);

    Vector<double> r0_icrs(3,0.0);
    Vector<double> v0_icrs(3,0.0);

    // (r, v, dr/dr0, dr/dv0, dv/dr0, dv/dv0, dr/dp0, dv/dp0)
    Vector<double> orbit0(42+6*numSRP,0.0);
    orbit0( 6) = 1.0; orbit0(10) = 1.0; orbit0(14) = 1.0;   // dr0/dr0
    orbit0(33) = 1.0; orbit0(37) = 1.0; orbit0(41) = 1.0;   // dv0/dv0

    satVectorMap satOrbit0;

    string block;
    map<SatID,string> satBlock;

    for(int i=1; i<=maxNumSat; ++i)
    {

//        if(i==1 || i==2 || i==4 || i==6 || i==11 || i==21) continue;

        if(i <= MAX_PRN_GPS)
        {
            sat.id = i;
            sat.system = SatID::systemGPS;
        }
        else if(i <= MAX_PRN_GPS+MAX_PRN_GAL)
        {
            sat.id = i - MAX_PRN_GPS;
            sat.system = SatID::systemGalileo;
        }
        else
        {
            sat.id = i - MAX_PRN_GPS - MAX_PRN_GAL;
            sat.system = SatID::systemBDS;
        }

        t2cRaw = refSys.T2CMatrix(utc0);
        t2cDot = refSys.dT2CMatrix(utc0);

        try
        {
            r0_itrs = sp3Store.getXvt(sat,gps0).x.toVector();
            v0_itrs = sp3Store.getXvt(sat,gps0).v.toVector();
        }
        catch(...)
        {
//            cerr << "No initial PV for " << sat << "." << endl;
            continue;
        }

        r0_icrs = t2cRaw * r0_itrs;
        v0_icrs = t2cRaw * v0_itrs + t2cDot * r0_itrs;

        // r0, v0
        orbit0(0) = r0_icrs(0); orbit0(1) = r0_icrs(1); orbit0(2) = r0_icrs(2);
        orbit0(3) = v0_icrs(0); orbit0(4) = v0_icrs(1); orbit0(5) = v0_icrs(2);

        satOrbit0[sat] = orbit0;
    }


    // Orbit Fit
    bool forward(false);

    cout << "Start Orbit Fit: " << endl;

    CommonTime gps(gps0);
    CommonTime utc(utc0);
    CommonTime tt(tt0);

    satVectorMap satOrbit(satOrbit0);

    vector< CommonTime > times;
    times.push_back(tt0);

    vector<satVectorMap> satOrbits;
    satOrbits.push_back(satOrbit0);


    double dt(0.0);

    Vector<double> orbit(42+6*numSRP,0.0);
    Vector<double> r_obs(3,0.0);


    epochMatrixMap epochOrbit;
    satEpochMatrixMap satEpochOrbit;

    Matrix<double> equation(3,8+numSRP,0.0);


    // rkf78
    rkf78.setCurrentTime(tt0);
    rkf78.setCurrentState(satOrbit0);

//    cout << fixed;

    for(int i=0; i<5*8; ++i)
    {
        if(forward) gps = gps + 60.0;
        else        gps = gps - 60.0;

        utc = refSys.GPS2UTC(gps);
        tt = refSys.GPS2TT(gps);

//        cout << CivilTime(gps) << endl;

        satOrbit = rkf78.integrateTo(tt);

        rkf78.setCurrentTime(tt);
        rkf78.setCurrentState(satOrbit);

        dt = gps - gps0;

        if( int(dt)%int(300.0) == 0 )
        {
            times.push_back(tt);
            satOrbits.push_back(satOrbit);
        }

        if( int(dt)%int(arcInt) != 0 ) continue;

        t2cRaw = refSys.T2CMatrix(utc);

        for(satVectorMap::const_iterator it = satOrbit.begin();
            it != satOrbit.end();
            ++it)
        {
            sat = it->first;
            orbit = it->second;

            try
            {
                r_obs = t2cRaw * sp3Store.getXvt(sat,gps).x.toVector();
            }
            catch(...)
            {
                cerr << "get sp3 orbit of " << sat << " error." << endl;
                continue;
            }

            for(int m=0; m<3; ++m)
            {
                equation(m,0) = r_obs(m);
                equation(m,1) = orbit(m);

                for(int n=0; n<3; ++n)
                {
                    equation(m,2+n) = orbit( 6+3*m+n);
                    equation(m,5+n) = orbit(15+3*m+n);
                }

                for(int n=0; n<numSRP; ++n)
                {
                    equation(m,8+n) = orbit(42+numSRP*m+n);
                }
            }

            satEpochMatrixMap::iterator it2( satEpochOrbit.find(sat) );

            // if sat exists
            if(it2 != satEpochOrbit.end())
            {
                epochOrbit = it2->second;

                epochMatrixMap::iterator it3( epochOrbit.find(gps) );

                // if epoch exists
                if(it3 != epochOrbit.end())
                {
                    continue;
                }
                // if epoch not exist, add it
                else
                {
                    epochOrbit[gps] = equation;
                    it2->second = epochOrbit;
                }
            }
            // if sat not exist
            else
            {
                epochOrbit[gps] = equation;
                satEpochOrbit[sat] = epochOrbit;
            }
        }

    } // End of 'for(int i=0;...)'


    // adams
    adams.setCurrentTime(times);
    adams.setCurrentState(satOrbits);

    vector<SatID> invalidSats;
    CommonTime gps_temp, utc_temp, tt_temp;
    satVectorMap invalidSatOrbit;

    while(true)
    {
        gps_temp = gps;
        utc_temp = refSys.GPS2UTC(gps_temp);
        tt_temp = refSys.GPS2TT(gps_temp);

        if(forward) gps = gps + 300.0;
        else        gps = gps - 300.0;

        utc = refSys.GPS2UTC(gps);
        tt = refSys.GPS2TT(gps);

        satOrbit = adams.integrateTo(tt);

        invalidSats = adams.getInvalidSats();

//        cout << CivilTime(gps) << ' ' << invalidSats.size() << endl;

        if(invalidSats.size() != 0)
        {
            // prepare initial state of invalid sats
            for(int i=0; i<invalidSats.size(); ++i)
            {
                sat = invalidSats[i];

                invalidSatOrbit[sat] = (satOrbits[8])[sat];
            }

            // set initial time and state of rkf78
            rkf78.setCurrentTime(tt_temp);
            rkf78.setCurrentState(invalidSatOrbit);

            // perform rkf78
            for(int j=0; j<5; ++j)
            {
                if(forward) gps_temp = gps_temp + 60.0;
                else        gps_temp = gps_temp - 60.0;

                utc_temp = refSys.GPS2UTC(gps_temp);
                tt_temp = refSys.GPS2TT(gps_temp);

                invalidSatOrbit = rkf78.integrateTo(tt_temp);

                rkf78.setCurrentTime(tt_temp);
                rkf78.setCurrentState(invalidSatOrbit);
            }

            // replace state of invalid sats
            for(int i=0; i<invalidSats.size(); ++i)
            {
                sat = invalidSats[i];

                satOrbit[sat] = invalidSatOrbit[sat];
            }
        }

        times.push_back(tt);
        times.erase(times.begin());
        satOrbits.push_back(satOrbit);
        satOrbits.erase(satOrbits.begin());

        adams.setCurrentTime(times);
        adams.setCurrentState(satOrbits);

        dt = gps - gps0;

        if(std::abs(dt) >= arcLen*3600.0) break;

        if( int(dt)%int(arcInt) != 0 ) continue;

        t2cRaw = refSys.T2CMatrix(utc);

        for(satVectorMap::const_iterator it = satOrbit.begin();
            it != satOrbit.end();
            ++it)
        {
            sat = it->first;
            orbit = it->second;

            try
            {
                r_obs = t2cRaw * sp3Store.getXvt(sat,gps).x.toVector();
            }
            catch(...)
            {
                cerr << "get sp3 orbit of " << sat << " error." << endl;
                continue;
            }

            for(int m=0; m<3; ++m)
            {
                equation(m,0) = r_obs(m);
                equation(m,1) = orbit(m);

                for(int n=0; n<3; ++n)
                {
                    equation(m,2+n) = orbit( 6+3*m+n);
                    equation(m,5+n) = orbit(15+3*m+n);
                }

                for(int n=0; n<numSRP; ++n)
                {
                    equation(m,8+n) = orbit(42+numSRP*m+n);
                }
            }

            satEpochMatrixMap::iterator it2( satEpochOrbit.find(sat) );

            // if sat exists
            if(it2 != satEpochOrbit.end())
            {
                epochOrbit = it2->second;

                epochMatrixMap::iterator it3( epochOrbit.find(gps) );

                // if epoch exists
                if(it3 != epochOrbit.end())
                {
                    continue;
                }
                // if epoch not exist, add it
                else
                {
                    epochOrbit[gps] = equation;
                    it2->second = epochOrbit;
                }
            }
            // if sat not exist
            else
            {
                epochOrbit[gps] = equation;
                satEpochOrbit[sat] = epochOrbit;
            }
        }

    } // End of 'while(true)'

//    cout << Counter::now() << endl;

    string icsFile;

    try
    {
        icsFile = confData.getValue("icsFileName", "DEFAULT");
    }
    catch(...)
    {
        cerr << "ICs File Name Get Error." << endl;
        return 1;
    }

    ofstream fics(icsFile.c_str());
    fics << fixed;

    for(satEpochMatrixMap::const_iterator it = satEpochOrbit.begin();
        it != satEpochOrbit.end();
        ++it)
    {
        sat = it->first;
        epochOrbit = it->second;

        int numEqu( epochOrbit.size() );

        Vector<double> omc(numEqu*3,0.0);
        Matrix<double> design(numEqu*3,6+numSRP,0.0);

        int count(0);
        for(epochMatrixMap::const_iterator it2 = epochOrbit.begin();
            it2 != epochOrbit.end();
            ++it2)
        {
            equation = it2->second;

            for(int i=0; i<3; ++i)
            {
                omc(3*count+i) = equation(i,0) - equation(i,1);

                for(int j=2; j<8+numSRP; ++j)
                {
                    design(3*count+i,j-2) = equation(i,j);
                }
            }

            ++count;
        }

//        cout << fixed << setprecision(6);
//        for(int i=0; i<3*numEqu; ++i)
//        {
//            cout << setw(20) << omc(i);
//
//            for(int j=0; j<6+numSRP; ++j)
//            {
//                cout << setw(20) << design(i,j);
//            }
//            cout << endl;
//        }
//
//        cout << endl;

//        continue;

        Matrix<double> A( transpose(design)*design );
        Vector<double> B( transpose(design)*omc );

        Vector<double> dx( inverse(A)*B );

        double sigma(0.0);

        omc -= design*dx;

        for(int i=1; i<=numEqu*3; ++i)
            sigma += omc(i-1)*omc(i-1);

        sigma = std::sqrt(sigma/(numEqu*3));

//        if(arcLen > 40 && sigma > 0.05) continue;
//        if(arcLen < 40 && sigma > 0.02) continue;

        cout << "Postfit RMS of "
             << sat << ": "
             << setw(10)
             << asString(sigma*1e2,3)
             << " (cm)" << endl;

        orbit0 = satOrbit0[sat];

        fics << sat;

        fics << fixed
             << setprecision(6)
             << setw(20) << dx( 0) + orbit0(0)
             << setw(20) << dx( 1) + orbit0(1)
             << setw(20) << dx( 2) + orbit0(2)
             << setprecision(6)
             << setw(15) << dx( 3) + orbit0(3)
             << setw(15) << dx( 4) + orbit0(4)
             << setw(15) << dx( 5) + orbit0(5);

        for(int i=0; i<numSRP; ++i)
        {
            fics << setprecision(3)
                 << setw(10) << dx(6+i) + srpc0(i);
        }

        fics << endl;
    }

//    cout << Counter::now() << endl;


    fics.close();

    return 0;

}
