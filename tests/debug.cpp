#include <iostream>
#include "Legendre.hpp"
#include "EGM08Model.hpp"
#include "constants.hpp"


using namespace std;
using namespace gpstk;


int main(void)
{
//    Vector<double> leg0, leg1, leg2;
//
//    double lat(30.0*DEG_TO_RAD);
//
//    legendre(12, lat, leg0, leg1, leg2, 2);
//
//    for(int i=0; i<leg0.size(); ++i)
//    {
//        cout << setw(15) << leg0(i) << endl;
//    }


    // EOP File
    EOPDataStore2 eopDataStore;
    try
    {
        string eopFile("/home/kfkuang/new/ROCKET/tables/finals2000A.data");
        eopDataStore.loadIERSFile(eopFile);
    }
    catch(...)
    {
        cerr << "EOP File Load Error." << endl;
        return 1;
    }

    eopDataStore.setInterpPoints(4);
    eopDataStore.setRegularization(false);


    // LeapSecond File
    LeapSecStore leapSecStore;
    try
    {
        string lsFile("/home/kfkuang/new/ROCKET/tables/Leap_Second.dat");
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
        string ephFile("/home/kfkuang/new/ROCKET/tables/1980_2040.DE405");
        solSys.initializeWithBinaryFile(ephFile);
    }
    catch(...)
    {
        cerr << "Solar System Initialize Error." << endl;
        return 1;
    }


    CivilTime time(2016,1,8,0,0,0.0, TimeSystem::UTC);
    CommonTime utc( time.convertToCommonTime() );
    CommonTime ut1( refSys.UTC2UT1(utc) );
    CommonTime tt( refSys.UTC2TT(utc) );

    Matrix<double> dCS;

    // solid tide
    EarthSolidTide solidTide;
    solidTide.setReferenceSystem(refSys);
    solidTide.setSolarSystem(solSys);

    dCS = solidTide.getSolidTide(tt);

//    cout << fixed << setprecision(15);
//    for(int i=0; i<dCS.rows(); ++i)
//    {
//        cout << setw(20) << dCS(i,0)
//             << setw(20) << dCS(i,1)
//             << endl;
//    }
//    cout << endl;

    // ocean tide
    EarthOceanTide oceanTide(8,8);
    oceanTide.setReferenceSystem(refSys);

    string otFileName("/home/kfkuang/new/ROCKET/tables/fes2004_Cnm-Snm.dat");

    try
    {
        oceanTide.loadFile(otFileName);
    }
    catch(...)
    {
        cerr << "ocean tide file '" << otFileName << "' load error." << endl;
        exit(-1);
    }

    dCS = oceanTide.getOceanTide(tt);

//    cout << fixed << setprecision(15);
//    for(int i=0; i<dCS.rows(); ++i)
//    {
//        cout << setw(20) << dCS(i,0)
//             << setw(20) << dCS(i,1)
//             << endl;
//    }
//    cout << endl;

    // pole tide
    EarthPoleTide poleTide;
    poleTide.setReferenceSystem(refSys);

    dCS = poleTide.getPoleTide(tt);

//    cout << fixed << setprecision(15);
//    for(int i=0; i<dCS.rows(); ++i)
//    {
//        cout << setw(20) << dCS(i,0)
//             << setw(20) << dCS(i,1)
//             << endl;
//    }
//    cout << endl;


    double BETA[6] = {0.0};
    double FNUT[5] = {0.0};

    DoodsonArguments(ut1,tt, BETA,FNUT);

    for(int i=0; i<6; ++i)
    {
        cout << setw(10) << BETA[i];
    }
    cout << endl;

    for(int i=0; i<5; ++i)
    {
        cout << setw(10) << FNUT[i];
    }
    cout << endl;

//    double jd = JulianDate(tt).jd;
//    SolarSystem::Planet center(SolarSystem::Earth);
//    SolarSystem::Planet target;
//
//    // moon position and velocity in ECI, unit: km, km/day
//    target = SolarSystem::Moon;
//    double rv_moon[6] = {0.0};
//    solSys.computeState(jd, target, center, rv_moon);
//
//    // moon position in ECI, unit: m
//    Vector<double> rm_eci(3,0.0);
//    rm_eci(0) = rv_moon[0];
//    rm_eci(1) = rv_moon[1];
//    rm_eci(2) = rv_moon[2];
//    rm_eci *= 1000.0;
//
//    cout << fixed << setprecision(2);
//
//    cout << setw(20) << rm_eci(0)
//         << setw(20) << rm_eci(1)
//         << setw(20) << rm_eci(2)
//         << endl;

    return 0;
}
