#pragma ident "$Id$"

#include "ConfDataReader.hpp"

#include "DataStructures.hpp"

#include "EOPDataStore2.hpp"

#include "LeapSecStore.hpp"

#include "ReferenceSystem.hpp"

#include "SP3EphemerisStore.hpp"

#include "Rinex3EphemerisStore2.hpp"

#include "SolarSystem.hpp"

#include "GNSSOrbit.hpp"

#include "RKF78Integrator.hpp"

#include "MSCStore.hpp"

#include "DCBDataReader.hpp"

#include "AntexReader.hpp"

#include "BLQDataReader.hpp"

#include "TropModel.hpp"

#include "CC2NONCC.hpp"

#include "RequireObservables.hpp"

#include "ConvertObservables.hpp"

#include "LinearCombinations.hpp"
#include "ComputeLinear.hpp"

#include "MWCSDetector.hpp"
#include "LICSDetector.hpp"

#include "SatArcMarker.hpp"
#include "SatArcMarker2.hpp"

#include "BasicModel.hpp"

#include "ComputeElevWeights.hpp"

#include "ComputeTropModel.hpp"

#include "ComputeStaTides.hpp"

#include "CorrectObservables.hpp"

#include "GravitationalDelay.hpp"

#include "ComputeSatAttitude.hpp"

#include "ComputeSatPCenter.hpp"

#include "ComputeWindUp.hpp"

#include "PhaseCodeAlignment.hpp"

#include "ComputeStaClock.hpp"

#include "ComputePartials.hpp"

#include "Epoch.hpp"

#include "Counter.hpp"

#include "Variable.hpp"

#include <Eigen/Dense>


using namespace std;
using namespace gpstk;
using namespace gpstk::StringUtils;


struct varData
{
    varData() : gps(CommonTime::BEGINNING_OF_TIME),
                state(0.0), covar(0.0)
    {};

    CommonTime gps;
    double state;
    double covar;
};


/// Returns 0 on success.
int main(int argc, char *argv[])
{

    double clock_start( Counter::now() );

    double dt(30.0);

    // conf file
    string confFileName("gps_orbit_clock.conf");

    ConfDataReader confReader;

    try
    {
        confReader.open(confFileName);
    }
    catch(...)
    {
        cerr << "conf file '" << confFileName
             << "' open error." << endl;
        exit(-1);
    }


    SourceIDSet allSourceSet;
    SatIDSet allSatSet;

    // station list
    string staListFileName;

    try
    {
        staListFileName = confReader.getValue("staListFileName", "DEFAULT");
    }
    catch(...)
    {
        cerr << "station list file name get error." << endl;
        exit(-1);
    }

    ifstream stationList;

    try
    {
        stationList.open(staListFileName.c_str(), ios::in);

        string line, station;

        while( !stationList.eof() && stationList.good() )
        {
            getline(stationList, line);

            if( stationList.eof() ) break;
            if( stationList.bad() ) break;

            if(line.size() != 0)
            {
                station = line.substr(0,4);

                SourceID source;

                transform(station.begin(),station.end(),station.begin(),::toupper);

                source = SourceID(SourceID::Mixed, station);

                allSourceSet.insert(source);
            }
        }
    }
    catch(...)
    {
        cerr << "station list file '" << staListFileName
             << "' open error." << endl;
        exit(-1);
    }


    // satellite list
    string satListFileName;

    try
    {
        satListFileName = confReader.getValue("satListFileName", "DEFAULT");
    }
    catch(...)
    {
        cerr << "satellite list file name get error." << endl;
        exit(-1);
    }

    ifstream satList;

    try
    {
        satList.open(satListFileName.c_str(), ios::in);

        string satStr;

        while( !satList.eof() && satList.good() )
        {
            getline(satList, satStr);

            if( satList.eof() || satList.bad() ) break;

            if(satStr.size() != 0)
            {
                SatID sat;

                string sys( satStr.substr(0,1) );
                int id( asInt(satStr.substr(1,2)) );

                if(sys == "G")
                {
                    sat = SatID(id, SatID::systemGPS);
                }
                else
                {
                    continue;
                }

                allSatSet.insert(sat);
            }
        }
    }
    catch(...)
    {
        cerr << "satellite list file '" << satListFileName
             << "' open error." << endl;
        exit(-1);
    }


    map<SourceID,Rinex3ObsStream*> sourceStreamMap;

    bool firstTime( true );

    bool updateEOP( false );
    bool updateLS( false );
    bool updateSP3( false );
    bool updateNAV( false );
    bool updateSNX( false );
    bool updateDCB( false );
    bool updateATX( false );
    bool updateBLQ( false );
    bool updateGPT2( false );
    bool updateOBS( false );


    // EOPDataStore2
    EOPDataStore2 eopStore;

    string eopFileName;

    try
    {
        eopFileName = confReader.getValue("eopFileName", "DEFAULT");
    }
    catch(...)
    {
        cerr << "eop file name get error." << endl;
        exit(-1);
    }

    try
    {
        eopStore.loadIERSFile( eopFileName );
    }
    catch(...)
    {
        cerr << "eop file '" << eopFileName << "' load error." << endl;
        exit(-1);
    }


    // LeapSecStore
    LeapSecStore lsStore;

    string lsFileName;

    try
    {
        lsFileName = confReader.getValue("lsFileName", "DEFAULT");
    }
    catch(...)
    {
        cerr << "ls file name get error." << endl;
        exit(-1);
    }

    try
    {
        lsStore.loadFile( lsFileName );
    }
    catch(...)
    {
        cerr << "ls file '" << lsFileName << "' load error." << endl;
        exit(-1);
    }


    // ReferenceSystem
    ReferenceSystem refSys;
    refSys.setEOPDataStore( eopStore );
    refSys.setLeapSecStore( lsStore );


    // SolarSystem
    SolarSystem solSys;

    string jplFileName;

    try
    {
        jplFileName = confReader.getValue("jplFileName", "DEFAULT");
    }
    catch(...)
    {
        cerr << "jpl file name get error." << endl;
        exit(-1);
    }

    try
    {
        solSys.initializeWithBinaryFile( jplFileName );
    }
    catch(...)
    {
        cerr << "solar system initialize error." << endl;
        exit(-1);
    }


    // EGM08Model
    EGM08Model egm(12,12);

    string egmFileName;

    try
    {
        egmFileName = confReader.getValue("egmFileName", "DEFAULT");
    }
    catch(...)
    {
        cerr << "egm file name get error." << endl;
        exit(-1);
    }

    try
    {
        egm.loadFile( egmFileName );
    }
    catch(...)
    {
        cerr << "egm file '" << egmFileName << "' load error." << endl;
        exit(-1);
    }

    egm.setReferenceSystem( refSys );


    // EarthSolidTide
    EarthSolidTide solidTide;
    solidTide.setReferenceSystem( refSys );
    solidTide.setSolarSystem( solSys );

    egm.setEarthSolidTide( solidTide );


    // EarthOceanTide
    EarthOceanTide oceanTide;
    oceanTide.setReferenceSystem( refSys );

    string otFileName;

    try
    {
        otFileName = confReader.getValue("otFileName", "DEFAULT");
    }
    catch(...)
    {
        cerr << "ot file name get error." << endl;
        exit(-1);
    }

    try
    {
        oceanTide.loadFile( otFileName );
    }
    catch(...)
    {
        cerr << "ocean tide file '" << otFileName << "' load error." << endl;
        exit(-1);
    }

    egm.setEarthOceanTide( oceanTide );


    // EarthPoleTide
    EarthPoleTide poleTide;
    poleTide.setReferenceSystem( refSys );

    egm.setEarthPoleTide( poleTide );


    // ThirdBody
    ThirdBody thd;
    thd.setReferenceSystem( refSys );
    thd.setSolarSystem( solSys );
    thd.enableAllPlanets();


    // ECOM1Model
    ECOM1Model srp;
    srp.setReferenceSystem( refSys );
    srp.setSolarSystem( solSys );

    int numSRPC(5);


    // Relativity
    Relativity rel;


    // GNSSOrbit
    GNSSOrbit gnss;
    gnss.setEGMModel( egm );
    gnss.setThirdBody( thd );
    gnss.setSRPModel( srp );
    gnss.setRelativity( rel );


    // RKF78Integrator
    RKF78Integrator rkf78;
    rkf78.setStepSize( dt );


    // SP3EphemerisStore
    SP3EphemerisStore sp3Store;


    // MSCStore
    MSCStore mscStore;


    // DCBDataReader
    DCBDataReader dcbReader;


    // AntexReader
    AntexReader antexReader;

    string atxFileName;

    try
    {
        atxFileName = confReader.getValue("atxFileName", "DEFAULT");
    }
    catch(...)
    {
        cerr << "get atx file name error." << endl;
        exit(-1);
    }

    try
    {
        antexReader.open( atxFileName );
    }
    catch(...)
    {
        cerr << "atx file '" << atxFileName << "' open error." << endl;
        exit(-1);
    }


    // BLQDataReader
    BLQDataReader blqReader;

    string blqFileName;

    try
    {
        blqFileName = confReader.getValue("blqFileName", "DEFAULT");
    }
    catch(...)
    {
        cerr << "get blq file name error." << endl;
        exit(-1);
    }

    try
    {
        blqReader.open( blqFileName );
    }
    catch(...)
    {
        cerr << "blq file '" << blqFileName << "' open error." << endl;
        exit(-1);
    }


    // ViennaTropModel
    ViennaTropModel viennaTM;

    string gpt2FileName;

    try
    {
        gpt2FileName = confReader.getValue("gpt2FileName", "DEFAULT");
    }
    catch(...)
    {
        cerr << "get gpt2 file name error." << endl;
        exit(-1);
    }

    try
    {
        viennaTM.loadFile( gpt2FileName );
    }
    catch(...)
    {
        cerr << "gpt2 file '" << gpt2FileName << "' load error." << endl;
        exit(-1);
    }


    map<SourceID, Position> sourcePositionMap;
    map<SourceID, Triple> sourceMonumentMap;
    map<SourceID, Antenna> sourceAntennaMap;


    // CC2NONCC
    CC2NONCC cc2noncc;


    // RequireObservables
    RequireObservables requireObs;
    requireObs.addRequiredType(SatID::systemGPS, TypeID::C1W);
    requireObs.addRequiredType(SatID::systemGPS, TypeID::C2W);
    requireObs.addRequiredType(SatID::systemGPS, TypeID::L1C);
    requireObs.addRequiredType(SatID::systemGPS, TypeID::L2W);


    // ConvertObservables
    ConvertObservables convertObs;


    // LinearCombinations
    LinearCombinations linearComb;


    // ComputeLinear, Cycle Slip
    ComputeLinear linearCS;
    linearCS.addLinear(SatID::systemGPS, linearComb.mw12CombOfGPS);
    linearCS.addLinear(SatID::systemGPS, linearComb.li12CombOfGPS);

    // ComputeLinear, Phase Code Alignment
    ComputeLinear linearAlign;
    linearAlign.addLinear(SatID::systemGPS, linearComb.q1CombOfGPSL1L2);
    linearAlign.addLinear(SatID::systemGPS, linearComb.q2CombOfGPSL1L2);


    // ComputeLinear, PC
    ComputeLinear linearPC;
    linearPC.addLinear(SatID::systemGPS, linearComb.pc12CombOfGPS);

    // ComputeLinear, LC
    ComputeLinear linearLC;
    linearLC.addLinear(SatID::systemGPS, linearComb.lc12CombOfGPS);


    // ComputeLinear, prefit PC
    ComputeLinear prefitPCForSPP;
    prefitPCForSPP.addLinear(SatID::systemGPS, linearComb.pc12PrefitOfGPS);


    // ComputeLinear, prefit PC for POD
    ComputeLinear prefitPCForPOD;
    prefitPCForPOD.addLinear(SatID::systemGPS, linearComb.pc12PrefitOfGPSForPOD);

    // ComputeLinear, prefit LC for POD
    ComputeLinear prefitLCForPOD;
    prefitLCForPOD.addLinear(SatID::systemGPS, linearComb.lc12PrefitOfGPSForPOD);


    // MWCSDetector
    MWCSDetector markCSMW;
    markCSMW.setCheckExist(true);
    markCSMW.setSatSystem(SatID::systemGPS);
    markCSMW.setObsType(TypeID::MW12);
    markCSMW.setLLIType1(TypeID::LLI1);
    markCSMW.setLLIType2(TypeID::LLI2);
    markCSMW.setResultType1(TypeID::CSL1);
    markCSMW.setResultType2(TypeID::CSL2);
    markCSMW.setDeltaTMax(61.0);
    markCSMW.setMaxNumLambdas(3.0*WL_WAVELENGTH_GPS_L1L2);
    markCSMW.setUseLLI(true);


    // LICSDetector
    LICSDetector markCSLI;
    markCSLI.setCheckExist(true);
    markCSLI.setSatSystem(SatID::systemGPS);
    markCSLI.setObsType(TypeID::LI12);
    markCSLI.setLLIType1(TypeID::LLI1);
    markCSLI.setLLIType2(TypeID::LLI2);
    markCSLI.setResultType1(TypeID::CSL1);
    markCSLI.setResultType2(TypeID::CSL2);
    markCSLI.setDeltaTMax(61.0);
    markCSLI.setMinThreshold(0.240);
    markCSLI.setLIDrift(0.002);
    markCSLI.setUseLLI(true);


    // SatArcMarker
    SatArcMarker markArc;
    markArc.setDeleteUnstableSats(true);
    markArc.setUnstablePeriod(30.0*5+1.0);


    // BasicModel
    BasicModel basicModel;
    basicModel.setMinElev(10.0);
    basicModel.setDefaultObs(SatID::systemGPS, TypeID::PC12);


    // ComputeElevWeights
    ComputeElevWeights elevWeights;


    // ComputeTropModel
    ComputeTropModel computeTM;
    computeTM.setTropModel(viennaTM);


    // ComputeStaTides
    ComputeStaTides staTides;
    staTides.setBLQDataReader(blqReader);
    staTides.setReferenceSystem(refSys);
    staTides.setSolarSystem(solSys);


    // CorrectObservables
    CorrectObservables correctObs;
    correctObs.setTideCorr(staTides);


    // GravitationalDelay
    GravitationalDelay gravDelay;


    // ComputeStaClock
    ComputeStaClock staClock;
    staClock.setMinSatOfGPS(4);
    staClock.setMaxRMSOfGPS(3.0);


    // ComputeSatAttitude
    ComputeSatAttitude satAttitude;
    satAttitude.setReferenceSystem(refSys);
    satAttitude.setSolarSystem(solSys);
    satAttitude.setAntexReader(antexReader);


    // ComputeSatPCenter
    ComputeSatPCenter satPCenter;
    satPCenter.setReferenceSystem(refSys);
    satPCenter.setSolarSystem(solSys);
    satPCenter.setAntexReader(antexReader);


    // ComputeWindUp
    ComputeWindUp windUp;
    windUp.setReferenceSystem(refSys);
    windUp.setSolarSystem(solSys);


    // PhaseCodeAlignment
    PhaseCodeAlignment phaseAlignL1;
    phaseAlignL1.setSatSystem(SatID::systemGPS);
    phaseAlignL1.setCodeType(TypeID::Q1OfL1L2);
    phaseAlignL1.setPhaseType(TypeID::L1);
    phaseAlignL1.setPhaseWavelength(L1_WAVELENGTH_GPS);
    phaseAlignL1.setUseSatArc(true);

    PhaseCodeAlignment phaseAlignL2;
    phaseAlignL2.setSatSystem(SatID::systemGPS);
    phaseAlignL2.setCodeType(TypeID::Q2OfL1L2);
    phaseAlignL2.setPhaseType(TypeID::L2);
    phaseAlignL2.setPhaseWavelength(L2_WAVELENGTH_GPS);
    phaseAlignL2.setUseSatArc(true);


    // ComputePartials
    ComputePartials partials;
    partials.setReferenceSystem(refSys);


    Matrix<double> t2cRaw(3,3,0.0), t2cDot(3,3,0.0);


    SourceIDSet sourceSet;
    SatIDSet satSet;

    int numSource(0), numSat(0);

    sourceVectorMap sourcePosECI;

    satVectorMap satOrbit;
    satVectorMap satSRPC;
    satValueMap satClock;


    satVectorMap oldSatSRPC;
    satValueMap oldSatClock;


    VariableVector prevVarVec, currVarVec;

    Eigen::VectorXd prevState, currState, tempState;
    Eigen::MatrixXd prevCovar, currCovar, tempCovar;


    // RX related variables
    Variable clkSta(TypeID::dcdtSta);
    clkSta.setSourceIndexed( true );
    clkSta.setSatIndexed( false );
    clkSta.setInitialVariance( 1e+2*1e+2 );

    Variable zwdSta(TypeID::wetMap);
    zwdSta.setSourceIndexed( true );
    zwdSta.setSatIndexed( false );
    zwdSta.setInitialVariance( 5e-1*5e-1 );

    VariableVector rxVarVec;
    rxVarVec.push_back( clkSta );
    rxVarVec.push_back( zwdSta );


    // SV related variables
    Variable rxSat(TypeID::dSatX);
    rxSat.setSourceIndexed( false );
    rxSat.setSatIndexed( true );
    rxSat.setInitialVariance( 1e+0*1e+0 );

    Variable rySat(TypeID::dSatY);
    rySat.setSourceIndexed( false );
    rySat.setSatIndexed( true );
    rySat.setInitialVariance( 1e+0*1e+0 );

    Variable rzSat(TypeID::dSatZ);
    rzSat.setSourceIndexed( false );
    rzSat.setSatIndexed( true );
    rzSat.setInitialVariance( 1e+0*1e+0 );

    Variable vxSat(TypeID::dSatVX);
    vxSat.setSourceIndexed( false );
    vxSat.setSatIndexed( true );
    vxSat.setInitialVariance( 1e-2*1e-2 );

    Variable vySat(TypeID::dSatVY);
    vySat.setSourceIndexed( false );
    vySat.setSatIndexed( true );
    vySat.setInitialVariance( 1e-2*1e-2 );

    Variable vzSat(TypeID::dSatVZ);
    vzSat.setSourceIndexed( false );
    vzSat.setSatIndexed( true );
    vzSat.setInitialVariance( 1e-2*1e-2 );

    Variable d0Sat(TypeID::dSRPC1);
    d0Sat.setSourceIndexed( false );
    d0Sat.setSatIndexed( true );
    d0Sat.setInitialVariance( 5e-1*5e-1 );

    Variable dcSat(TypeID::dSRPC2);
    dcSat.setSourceIndexed( false );
    dcSat.setSatIndexed( true );
    dcSat.setInitialVariance( 5e-1*5e-1 );

    Variable dsSat(TypeID::dSRPC3);
    dsSat.setSourceIndexed( false );
    dsSat.setSatIndexed( true );
    dsSat.setInitialVariance( 5e-1*5e-1 );

    Variable y0Sat(TypeID::dSRPC4);
    y0Sat.setSourceIndexed( false );
    y0Sat.setSatIndexed( true );
    y0Sat.setInitialVariance( 5e-1*5e-1 );

    Variable ycSat(TypeID::dSRPC5);
    ycSat.setSourceIndexed( false );
    ycSat.setSatIndexed( true );
    ycSat.setInitialVariance( 5e-1*5e-1 );

    Variable ysSat(TypeID::dSRPC6);
    ysSat.setSourceIndexed( false );
    ysSat.setSatIndexed( true );
    ysSat.setInitialVariance( 5e-1*5e-1 );

    Variable b0Sat(TypeID::dSRPC7);
    b0Sat.setSourceIndexed( false );
    b0Sat.setSatIndexed( true );
    b0Sat.setInitialVariance( 5e-1*5e-1 );

    Variable bcSat(TypeID::dSRPC8);
    bcSat.setSourceIndexed( false );
    bcSat.setSatIndexed( true );
    bcSat.setInitialVariance( 5e-1*5e-1 );

    Variable bsSat(TypeID::dSRPC9);
    bsSat.setSourceIndexed( false );
    bsSat.setSatIndexed( true );
    bsSat.setInitialVariance( 5e-1*5e-1 );

    Variable clkSat(TypeID::dcdtSat);
    clkSat.setSourceIndexed( false );
    clkSat.setSatIndexed( true );
    clkSat.setInitialVariance( 1e+2*1e+2 );

    VariableVector svVarVec;
    svVarVec.push_back( rxSat );
    svVarVec.push_back( rySat );
    svVarVec.push_back( rzSat );
    svVarVec.push_back( vxSat );
    svVarVec.push_back( vySat );
    svVarVec.push_back( vzSat );
    svVarVec.push_back( d0Sat );
    svVarVec.push_back( dcSat );
    svVarVec.push_back( dsSat );
    svVarVec.push_back( y0Sat );
//    svVarVec.push_back( ycSat );
//    svVarVec.push_back( ysSat );
    svVarVec.push_back( b0Sat );
//    svVarVec.push_back( bcSat );
//    svVarVec.push_back( bsSat );
    svVarVec.push_back( clkSat );

    // RX-SV related variables
    Variable ambiVar(TypeID::BLC);
    ambiVar.setSourceIndexed( true );
    ambiVar.setSatIndexed( true );
    ambiVar.setInitialVariance( 1e+2*1e+2 );


    map<SourceID,int> sourceIndexMap;
    map<SatID,int> satIndexMap;
    map< SourceID, map<SatID,int> > sourceSatIndexMap;


    TypeIDSet keepTypes;

    keepTypes.insert(TypeID::prefitC12ForPOD);
    keepTypes.insert(TypeID::prefitL12ForPOD);
    keepTypes.insert(TypeID::weight);
    keepTypes.insert(TypeID::wetMap);
    keepTypes.insert(TypeID::CSL1);
    keepTypes.insert(TypeID::rho);



    // begin/break time
    string time1, time2;
    try
    {
        time1 = confReader.getValue("beginTime", "DEFAULT");
        time2 = confReader.getValue("breakTime", "DEFAULT");
    }
    catch(...)
    {
        cerr << "begin/break time get error." << endl;
        exit(-1);
    }

    int year, mon, day, hour, min;
    double sec;

    year = asInt( time1.substr( 0,4) );
    mon  = asInt( time1.substr( 5,2) );
    day  = asInt( time1.substr( 8,2) );
    hour = asInt( time1.substr(11,2) );
    min  = asInt( time1.substr(14,2) );
    sec  = asDouble( time1.substr(17,2) );

    CivilTime civil1(year,mon,day,hour,min,sec, TimeSystem::GPS);
    CommonTime gps1( civil1.convertToCommonTime() );

    cerr << "begin time: " << civil1 << endl;

    year = asInt( time2.substr( 0,4) );
    mon  = asInt( time2.substr( 5,2) );
    day  = asInt( time2.substr( 8,2) );
    hour = asInt( time2.substr(11,2) );
    min  = asInt( time2.substr(14,2) );
    sec  = asDouble( time2.substr(17,2) );

    CivilTime civil2(year,mon,day,hour,min,sec, TimeSystem::GPS);
    CommonTime gps2( civil2.convertToCommonTime() );

    cerr << "break time: " << civil2 << endl;


    // sp3 file dir
    string sp3FileDir;

    try
    {
        sp3FileDir = confReader.getValue("sp3FileDir", "DEFAULT");
    }
    catch(...)
    {
        cerr << "sp3 file directory get error." << endl;
        exit(-1);
    }

    // clk file dir
    string clkFileDir;

    try
    {
        clkFileDir = confReader.getValue("clkFileDir", "DEFAULT");
    }
    catch(...)
    {
        cerr << "clk file directory get error." << endl;
        exit(-1);
    }

    // analysis center
//    string analysisCenter;
//
//    try
//    {
//        analysisCenter = confReader.getValue("analysisCenter", "DEFAULT");
//    }
//    catch(...)
//    {
//        cerr << "analysis center get error." << endl;
//        exit(-1);
//    }

    // product abbr
    string productAbbr;

    try
    {
        productAbbr = confReader.getValue("productAbbr", "DEFAULT");
    }
    catch(...)
    {
        cerr << "product abbreviation get error." << endl;
        exit(-1);
    }

    // nav file dir
    string navFileDir;

    try
    {
        navFileDir = confReader.getValue("navFileDir", "DEFAULT");
    }
    catch(...)
    {
        cerr << "nav file directory get error." << endl;
        exit(-1);
    }

    // msc file dir
    string mscFileDir;

    try
    {
        mscFileDir = confReader.getValue("mscFileDir", "DEFAULT");
    }
    catch(...)
    {
        cerr << "msc file directory get error." << endl;
        exit(-1);
    }

    // dcb file dir
    string dcbFileDir;

    try
    {
        dcbFileDir = confReader.getValue("dcbFileDir", "DEFAULT");
    }
    catch(...)
    {
        cerr << "dcb file directory get error." << endl;
        exit(-1);
    }

    // obs file dir
    string obsFileDir;

    try
    {
        obsFileDir = confReader.getValue("obsFileDir", "DEFAULT");
    }
    catch(...)
    {
        cerr << "obs file directory get error." << endl;
        exit(-1);
    }

    // srpc noise
    Vector<double> srpcNoise(numSRPC,0.0);

    try
    {
        srpcNoise(0) = confReader.getValueAsDouble("d0Noise", "DEFAULT");
        srpcNoise(1) = confReader.getValueAsDouble("y0Noise", "DEFAULT");
        srpcNoise(2) = confReader.getValueAsDouble("b0Noise", "DEFAULT");
        srpcNoise(3) = confReader.getValueAsDouble("bcNoise", "DEFAULT");
        srpcNoise(4) = confReader.getValueAsDouble("bsNoise", "DEFAULT");
    }
    catch(...)
    {
        cerr << "srpc noise get error." << endl;
        exit(-1);
    }

    // num threads
    int numThreads;

    try
    {
        numThreads = confReader.getValueAsInt("numThreads", "DEFAULT");
    }
    catch(...)
    {
        cerr << "num of threads get error." << endl;
        exit(-1);
    }


    CommonTime gps( gps1 ), utc, tt;


    // orbit file name
    string orbitFileName;

    try
    {
        orbitFileName = confReader.getValue("orbitFileName", "DEFAULT");
    }
    catch(...)
    {
        cerr << "orbit file name get error." << endl;
        exit(-1);
    }

    // clock file name
    string clockFileName;

    try
    {
        clockFileName = confReader.getValue("clockFileName", "DEFAULT");
    }
    catch(...)
    {
        cerr << "clock file name get error." << endl;
        exit(-1);
    }


    // debug file name
    string debugFileName;

    try
    {
        debugFileName = confReader.getValue("debugFileName", "DEFAULT");
    }
    catch(...)
    {
        cerr << "debug file name get error." << endl;
        exit(-1);
    }


    cout << fixed;

    ofstream debugStream( debugFileName.c_str() );

    debugStream << fixed;

    map<SourceID, varData> sourceZWDData;

    while( true )
    {
        YDSTime yds( gps );
        GPSWeekSecond gws( gps );
        CivilTime ymd( gps );

        int doy( yds.doy );

        int year( ymd.year );
        int month( ymd.month );
        int day( ymd.day );
        int hour( ymd.hour );
        int minute( ymd.minute );
        double second( ymd.second );

        int week( gws.week );
        int dow( int(gws.sow/86400.0) );

        // if first epoch, update all files
        if( firstTime )
        {
            updateSP3 = true;
            updateNAV = true;
            updateSNX = true;
            updateDCB = true;
            updateOBS = true;

            firstTime = false;
        }


        // update SP3, every day
        if(hour==0 && minute==0 && second==0.0) updateSP3 = true;

        if( updateSP3 )
        {

            sp3Store = SP3EphemerisStore();

            sp3Store.rejectBadPositions( true );
            sp3Store.rejectBadClocks( true );
            sp3Store.setPosGapInterval( 900.0 + 1.0 );
            sp3Store.setPosMaxInterval( 9*900.0 + 1.0 );

            string sp3FileHdr(sp3FileDir + "igs");
            string clkFileHdr(clkFileDir + "igs");

            string sp3File, clkFile;

            try
            {
                // current day - 1
                if(dow == 0)
                {
                    sp3File = sp3FileHdr + asString(week-1) + "6" + ".sp3";
                    clkFile = clkFileHdr + asString(week-1) + "6" + ".clk_30s";
                }
                else
                {
                    sp3File = sp3FileHdr + asString(week) + asString(dow-1) + ".sp3";
                    clkFile = clkFileHdr + asString(week) + asString(dow-1) + ".clk_30s";
                }

                sp3Store.loadSP3File( sp3File );
                sp3Store.loadRinexClockFile( clkFile );

                // current day
                sp3File = sp3FileHdr + asString(week) + asString(dow) + ".sp3";
                clkFile = clkFileHdr + asString(week) + asString(dow) + ".clk_30s";

                sp3Store.loadSP3File( sp3File );
                sp3Store.loadRinexClockFile( clkFile );


                // current day + 1
                if(dow == 6)
                {
                    sp3File = sp3FileHdr + asString(week+1) + "0" + ".sp3";
                    clkFile = clkFileHdr + asString(week+1) + "0" + ".clk_30s";
                }
                else
                {
                    sp3File = sp3FileHdr + asString(week) + asString(dow+1) + ".sp3";
                    clkFile = clkFileHdr + asString(week) + asString(dow+1) + ".clk_30s";
                }

                sp3Store.loadSP3File( sp3File );
                sp3Store.loadRinexClockFile( clkFile );

            }
            catch(...)
            {
                break;
            }

            basicModel.setEphStore( sp3Store );
            satPCenter.setEphStore( sp3Store );
            windUp.setEphStore( sp3Store );

            updateSP3 = false;
        }


        // update SNX, every day
        if(hour==0 && minute==0 && second==0.0) updateSNX = true;

        if( updateSNX )
        {
            string file( mscFileDir );

            string yr2;

            if(month == 1 && day == 1)
            {
                yr2 = asString(year-1).substr(2,2);
            }
            else
            {
                yr2 = asString(year).substr(2,2);
            }

            if(dow == 0)
            {
                file += "igs" + yr2 + "P" + asString(week-1) + "6" + ".msc";
            }
            else
            {
                file += "igs" + yr2 + "P" + asString(week) + asString(dow-1) + ".msc";
            }

            try
            {
                mscStore.loadFile( file );
            }
            catch(...)
            {
                break;
            }

            basicModel.setMSCStore( mscStore );
            computeTM.setMSCStore( mscStore );
            correctObs.setMSCStore( mscStore );
            gravDelay.setMSCStore( mscStore );
            satPCenter.setMSCStore( mscStore );
            windUp.setMSCStore( mscStore );
            partials.setMSCStore( mscStore );

            updateSNX = false;
        }


        // update DCB, every month
        if(day==1 && hour==0 && minute==0 && second==0.0) updateDCB = true;

        if( updateDCB )
        {
            string file( dcbFileDir );

            string yr2( asString(year).substr(2,2) );

            if(month < 10)
            {
                file += "P1C1" + yr2 + "0" + asString(month) + ".DCB";
            }
            else
            {
                file += "P1C1" + yr2 + asString(month) + ".DCB";
            }

            try
            {
                dcbReader.open( file );
            }
            catch(...)
            {
                break;
            }

            cc2noncc.setDCBDataReader( dcbReader );

            updateDCB = false;
        }


        // obs data
        gnssDataMap gData;

        // update OBS, every day
        if(hour==0 && minute==0 && second==0.0) updateOBS = true;

        if( updateOBS )
        {
            sourceStreamMap.clear();

            sourcePositionMap.clear();
            sourceMonumentMap.clear();
            sourceAntennaMap.clear();

            RinexObsID GC1C("GC1C"), GC1W("GC1W"), GL1C("GL1C");
            RinexObsID GC2W("GC2W"), GL2W("GL2W");

            for(SourceIDSet::iterator it = allSourceSet.begin();
                it != allSourceSet.end();
                ++it)
            {
                SourceID source( *it );

                string station( source.sourceName.substr(0,4) );

                string lower( station ), upper( station );
                transform(lower.begin(),lower.end(),lower.begin(),::tolower);
                transform(upper.begin(),upper.end(),upper.begin(),::toupper);


                // check if station exists in SNX

                MSCData mscData;

                try
                {
                    CommonTime temp( gps );

                    temp.setTimeSystem( TimeSystem::Unknown );

                    mscData = mscStore.findMSC(upper,temp);
                }
                catch(...)
                {
                    cerr << CivilTime(gps);
                    cerr << " station " << upper << " not exist in SNX file." << endl;
                    continue;
                }


                Position position( mscData.coordinates );


                // check if station exists in BLQ
                if( !blqReader.isValid(upper) )
                {
                    cerr << CivilTime(gps);
                    cerr << " station " << upper << " not exist in BLQ file." << endl;
                    continue;
                }


                // obs file name
                string file( obsFileDir );

                yds = YDSTime( gps );

                string yr4( asString(year) );
                string yr2( yr4.substr(2,2) );
                string dy3;

                if(doy < 10)        dy3 = "00" + asString(doy);
                else if(doy < 100)  dy3 = "0" + asString(doy);
                else                dy3 = asString(doy);

                file += yr4 + dy3 + "/" + lower + dy3 + "0." + yr2 + "o";

                // obs stream
                Rinex3ObsStream* pObsStream;

                try
                {
                    pObsStream = new Rinex3ObsStream();

                    if( !pObsStream )
                    {
                        delete pObsStream;
                        continue;
                    }

                    pObsStream->exceptions(ios::failbit);
                    pObsStream->open(file, ios::in);

                    // obs header
                    Rinex3ObsHeader roh;

                    try
                    {
                        (*pObsStream) >> roh;
                    }
                    catch(...)
                    {
                        (*pObsStream).close(); continue;
                    }


                    // check if types exist
                    map< string, vector<RinexObsID> > mapObsTypes;
                    mapObsTypes = roh.mapObsTypes;

                    vector<RinexObsID> roi;

                    if(mapObsTypes.find("G") != mapObsTypes.end())
                    {
                        roi = mapObsTypes["G"];

                        if( find(roi.begin(),roi.end(), GC1C) == roi.end() &&
                            find(roi.begin(),roi.end(), GC1W) == roi.end() )
                        {
                            (*pObsStream).close(); continue;
                        }

                        if( find(roi.begin(),roi.end(), GC2W) == roi.end() ||
                            find(roi.begin(),roi.end(), GL1C) == roi.end() ||
                            find(roi.begin(),roi.end(), GL2W) == roi.end() )
                        {
                            (*pObsStream).close(); continue;
                        }
                    }

                    Triple monument( roh.antennaDeltaHEN );

                    // check if antenna exists
                    string antType( roh.antType );


                    Antenna antenna;

                    try
                    {
                        antenna = antexReader.getAntenna( antType );
                    }
                    catch(ObjectNotFound& notFound)
                    {
                        antType.replace(16,4,"NONE");
                        antenna = antexReader.getAntenna( antType );
                    }
                    catch(...)
                    {
                        (*pObsStream).close(); continue;
                    }

                    sourcePositionMap[source] = position;
                    sourceMonumentMap[source] = monument;
                    sourceAntennaMap[source] = antenna;

                }
                catch(...)
                {
                    delete pObsStream;
                    pObsStream = (Rinex3ObsStream*)0;

                    continue;
                }


                gnssRinex gRin;

                streampos sp( pObsStream->tellg() );

                try
                {
                    (*pObsStream) >> gRin;
                }
                catch(...)
                {
                    (*pObsStream).close(); continue;
                }

                if(gRin.header.epoch == gps)
                {
                    gData.addGnssRinex( gRin );
                }
                else if(gRin.header.epoch > gps)
                {
                    pObsStream->seekg( sp );
                }
                else if(gRin.header.epoch < gps)
                {
                    while( true )
                    {
                        if( (*pObsStream).eof() )
                        {
                            break;
                        }

                        if(gRin.header.epoch == gps)
                        {
                            gData.addGnssRinex( gRin );
                            break;
                        }

                        (*pObsStream) >> gRin;
                    }
                }

                sourceStreamMap[source] = pObsStream;

            } // End of 'for(...)'

            correctObs.setSourceMonument( sourceMonumentMap );
            correctObs.setSourceAntenna( sourceAntennaMap );

            updateOBS = false;
        }
        else
        {
            SourceIDSet sourceRejectedSet;

            for(map<SourceID,Rinex3ObsStream*>::iterator it = sourceStreamMap.begin();
                it != sourceStreamMap.end();
                ++it)
            {
                SourceID source( it->first );

                Rinex3ObsStream* pObsStream( it->second );

                gnssRinex gRin;

                streampos sp( pObsStream->tellg() );

                try
                {
                    (*pObsStream) >> gRin;
                }
                catch(...)
                {
                    sourceRejectedSet.insert( source );
                    (*pObsStream).close();
                    continue;
                }

                if(gRin.header.epoch == gps)
                {
                    gData.addGnssRinex( gRin );
                }
                else if(gRin.header.epoch > gps)
                {
                    pObsStream->seekg( sp );
                }
                else if(gRin.header.epoch < gps)
                {
                    while( true )
                    {
                        if( (*pObsStream).eof() )
                        {
                            sourceRejectedSet.insert( source );
                            break;
                        }

                        if(gRin.header.epoch == gps)
                        {
                            gData.addGnssRinex( gRin );
                            break;
                        }

                        (*pObsStream) >> gRin;
                    }
                }

            } // End of 'for(...)'

            for(SourceIDSet::const_iterator it = sourceRejectedSet.begin();
                it != sourceRejectedSet.end();
                ++it)
            {
                sourceStreamMap.erase( *it );
            }
        }


        if(gData.size() == 0) break;


        gData.keepOnlySatID( allSatSet );


        try
        {
            double clock1( Counter::now() );

            utc = refSys.GPS2UTC( gps );
            tt = refSys.GPS2TT( gps );

            t2cRaw = refSys.T2CMatrix(utc);
            t2cDot = refSys.dT2CMatrix(utc);

            SatIDSet satRejectedSet;

            // final sat orbit
            satVectorMap finalSatOrbitECI;
            satValueMap finalSatClock;

            for(SatIDSet::iterator it = allSatSet.begin();
                it != allSatSet.end();
                ++it)
            {
                SatID sat( *it );

                Vector<double> rsat_t(3,0.0), vsat_t(3,0.0);
                Vector<double> rsat_c(3,0.0), vsat_c(3,0.0);
                Vector<double> orbit(6,0.0);

                double clock(0.0);

                try
                {
                    rsat_t = sp3Store.getXvt(sat,gps).x.toVector();
                    vsat_t = sp3Store.getXvt(sat,gps).v.toVector();
                    clock = sp3Store.getXvt(sat,gps).clkbias*C_MPS;

                    rsat_c = t2cRaw * rsat_t;
                    vsat_c = t2cRaw * vsat_t + t2cDot * rsat_t;

                    for(int i=0; i<3; ++i)
                    {
                        orbit(i+0) = rsat_c(i);
                        orbit(i+3) = vsat_c(i);
                    }

                    finalSatOrbitECI[sat] = orbit;
                    finalSatClock[sat] = clock;
                }
                catch(...)
                {
                    satRejectedSet.insert(sat);
                    continue;
                }
            }

            gData.removeSatID(satRejectedSet);


            // preprocessing
            gData >> cc2noncc               // C1C/C2W --> C1W/C2W
                  >> requireObs             // C1W/C2W/L1C/L2W
                  >> convertObs             // C1W/C2W/L1C/L2W --> C1/C2/L1/L2

                  >> linearPC               // PC
                  >> basicModel             // r,v,cdt of SV
                                            // rho,elev,azim of RX-SV
                  >> elevWeights            // weight
                  >> computeTM              // troposphere
                  >> correctObs             // disp,ARP,PCO,PCV of RX
                  >> gravDelay              // gravitational delay
                  >> linearPC
                  >> prefitPCForSPP
                  >> staClock;

            sourceValueMap sourceClock;
            sourceClock = staClock.getSourceClockMap();

            partials.setSourceClockMap( sourceClock );

            gData >> linearCS               // Cycle Slip
                  >> markCSMW               // MW detector
                  >> markCSLI               // LI detector
                  >> markArc                // obs arc

                  >> partials;              // partials

            satYawDataMap satYawData;
            satYawData = satAttitude.getAttitude(gps,finalSatOrbitECI);

            satPCenter.setSatYawData( satYawData );
            windUp.setSatYawData( satYawData );

            gData >> satPCenter             // PCO,PCV of SV
                  >> windUp;                // Wind Up

            gData >> linearAlign            // Phase Code Align
                  >> phaseAlignL1           // Align L1
                  >> phaseAlignL2           // Align L2

                  >> linearPC               // PC
                  >> linearLC               // LC
                  >> prefitPCForPOD
                  >> prefitLCForPOD;

            sourcePosECI = partials.getSourcePosECIMap();

            double clock2( Counter::now() );

            gData.keepOnlyTypeID( keepTypes );


            map<SatID,int> satSourceNumMap;

            int numObs(0), numUnk(0);

            int iteration(0);

            while(iteration < 2)
            {
                iteration++;

                currVarVec.clear();

                sourceIndexMap.clear();
                satIndexMap.clear();
                sourceSatIndexMap.clear();


                //////// start of filter initialization ////////

                // number of source
                sourceSet = gData.getSourceIDSet();
                numSource = sourceSet.size();

                // number of sat
                satSet = gData.getSatIDSet();
                numSat = satSet.size();

                if(numSource == 0 || numSat == 0) break;

                satSourceNumMap.clear();

                for(gnssDataMap::iterator gdmIt = gData.begin();
                    gdmIt != gData.end();
                    ++gdmIt)
                {
                    for(sourceDataMap::iterator sdmIt = gdmIt->second.begin();
                        sdmIt != gdmIt->second.end();
                        ++sdmIt)
                    {
                        for(satTypeValueMap::iterator stvmIt = sdmIt->second.begin();
                            stvmIt != sdmIt->second.end();
                            ++stvmIt)
                        {
                            SatID sat( stvmIt->first );

                            if(satSourceNumMap.find(sat) != satSourceNumMap.end())
                            {
                                satSourceNumMap[sat] += 1;
                            }
                            else
                            {
                                satSourceNumMap[sat] = 1;
                            }
                        }
                    }
                }


                int id(0);

                // source-related variable
                for(SourceIDSet::iterator it = sourceSet.begin();
                    it != sourceSet.end();
                    ++it)
                {
                    SourceID source(*it);

                    sourceIndexMap[source] = id;

                    for(VariableVector::iterator varIt = rxVarVec.begin();
                        varIt != rxVarVec.end();
                        ++varIt)
                    {
                        Variable var(*varIt);
                        var.setSource( source );
                        var.setCurrentIndex( id );

                        currVarVec.push_back( var );

                        id++;
                    }
                }

                // sat-related variable
                for(SatIDSet::iterator it = satSet.begin();
                    it != satSet.end();
                    ++it)
                {
                    SatID sat(*it);

                    satIndexMap[sat] = id;

                    for(VariableVector::iterator varIt = svVarVec.begin();
                        varIt != svVarVec.end();
                        ++varIt)
                    {
                        Variable var(*varIt);
                        var.setSatellite( sat );
                        var.setCurrentIndex( id );

                        currVarVec.push_back( var );

                        id++;
                    }
                }

                // source-sat-related variable
                numObs = 0;
                for(gnssDataMap::iterator gdmIt = gData.begin();
                    gdmIt != gData.end();
                    ++gdmIt)
                {
                    for(sourceDataMap::iterator sdmIt = gdmIt->second.begin();
                        sdmIt != gdmIt->second.end();
                        ++sdmIt)
                    {
                        SourceID source( sdmIt->first );

                        for(satTypeValueMap::iterator stvmIt = sdmIt->second.begin();
                            stvmIt != sdmIt->second.end();
                            ++stvmIt)
                        {
                            SatID sat(stvmIt->first);
                            typeValueMap tvm(stvmIt->second);

                            if( tvm.find(TypeID::prefitC12ForPOD) == tvm.end() ||
                                tvm.find(TypeID::prefitL12ForPOD) == tvm.end() )
                                continue;

                            Variable var(ambiVar);
                            var.setSource( source );
                            var.setSatellite( sat );
                            var.setCurrentIndex( id );

                            currVarVec.push_back( var );

                            sourceSatIndexMap[source][sat] = id;

                            numObs += 2;
                            id++;
                        }
                    }
                }

                // set up index
                for(VariableVector::iterator currIt = currVarVec.begin();
                    currIt != currVarVec.end();
                    ++currIt)
                {
                    TypeID type( currIt->getType() );

                    VariableVector::iterator prevIt;
                    prevIt = find(prevVarVec.begin(),prevVarVec.end(), *currIt);

                    if(prevIt != prevVarVec.end())
                    {
                        currIt->setPreviousIndex( prevIt->getCurrentIndex() );

                        if(type == TypeID::BLC)
                        {
                            SourceID source( currIt->getSource() );
                            SatID sat( currIt->getSatellite() );

                            if(gData.getValue(gps,source,sat,TypeID::CSL1) != 0.0)
                            {
                                currIt->setPreviousIndex( -1 );
                            }
                        }
                    }
                }


                // state and covar
                numUnk = currVarVec.size();

                tempState.resize(numUnk);
                for(int i=0; i<numUnk; ++i) tempState(i) = 0.0;
                tempCovar = Eigen::MatrixXd::Zero(numUnk,numUnk);

                for(int i=0; i<numUnk; ++i)
                {
                    Variable var1( currVarVec[i] );

                    int prevIndex1( var1.getPreviousIndex() );
                    int currIndex1( var1.getCurrentIndex() );

                    TypeID type( var1.getType() );

                    // new variable
                    if(prevIndex1 == -1)
                    {
                        if(type == TypeID::dSatX)
                        {
                            SatID sat( var1.getSatellite() );

                            if(finalSatOrbitECI.find(sat) != finalSatOrbitECI.end())
                            {
                                tempState(currIndex1) = finalSatOrbitECI[sat](0);
                            }
                            tempCovar(currIndex1,currIndex1) = var1.getInitialVariance();
                        }
                        else if(type == TypeID::dSatY)
                        {
                            SatID sat( var1.getSatellite() );

                            if(finalSatOrbitECI.find(sat) != finalSatOrbitECI.end())
                            {
                                tempState(currIndex1) = finalSatOrbitECI[sat](1);
                            }
                            tempCovar(currIndex1,currIndex1) = var1.getInitialVariance();
                        }
                        else if(type == TypeID::dSatZ)
                        {
                            SatID sat( var1.getSatellite() );

                            if(finalSatOrbitECI.find(sat) != finalSatOrbitECI.end())
                            {
                                tempState(currIndex1) = finalSatOrbitECI[sat](2);
                            }
                            tempCovar(currIndex1,currIndex1) = var1.getInitialVariance();
                        }
                        else if(type == TypeID::dSatVX)
                        {
                            SatID sat( var1.getSatellite() );

                            if(finalSatOrbitECI.find(sat) != finalSatOrbitECI.end())
                            {
                                tempState(currIndex1) = finalSatOrbitECI[sat](3);
                            }
                            tempCovar(currIndex1,currIndex1) = var1.getInitialVariance();
                        }
                        else if(type == TypeID::dSatVY)
                        {
                            SatID sat( var1.getSatellite() );

                            if(finalSatOrbitECI.find(sat) != finalSatOrbitECI.end())
                            {
                                tempState(currIndex1) = finalSatOrbitECI[sat](4);
                            }
                            tempCovar(currIndex1,currIndex1) = var1.getInitialVariance();
                        }
                        else if(type == TypeID::dSatVZ)
                        {
                            SatID sat( var1.getSatellite() );

                            if(finalSatOrbitECI.find(sat) != finalSatOrbitECI.end())
                            {
                                tempState(currIndex1) = finalSatOrbitECI[sat](5);
                            }
                            tempCovar(currIndex1,currIndex1) = var1.getInitialVariance();
                        }
                        else if(type == TypeID::dcdtSat)
                        {
                            SatID sat( var1.getSatellite() );

                            if(finalSatClock.find(sat) != finalSatClock.end())
                            {
                                tempState(currIndex1) = finalSatClock[sat];
                                tempCovar(currIndex1,currIndex1) = 1e+2*1e+2;
                            }
                            else
                            {
                                tempState(currIndex1) = 0.0;
                                tempCovar(currIndex1,currIndex1) = 3e+5*3e+5;
                            }
                        }
                        else if(type == TypeID::dSRPC1)
                        {
                            SatID sat( var1.getSatellite() );

                            tempState(currIndex1) = -1.0;
                            tempCovar(currIndex1,currIndex1) = var1.getInitialVariance();
                        }
                        else if(type == TypeID::wetMap)
                        {
                            SourceID source( var1.getSource() );

                            if(sourceZWDData.find(source) != sourceZWDData.end())
                            {
                                CommonTime gps0(CommonTime::BEGINNING_OF_TIME);

                                double state(0.0), covar(0.0);

                                gps0 = sourceZWDData[source].gps;

                                state = sourceZWDData[source].state;
                                covar = sourceZWDData[source].covar;

                                covar += 3e-8*(gps-gps0);

                                if(covar > 5e-1*5e-1)
                                {
                                    state = 0.0;
                                    covar = 5e-1*5e-1;
                                }

                                tempState(currIndex1) = state;
                                tempCovar(currIndex1,currIndex1) = covar;
                            }
                            else
                            {
                                tempState(currIndex1) = 0.0;
                                tempCovar(currIndex1,currIndex1) = var1.getInitialVariance();
                            }
                        }
                        else
                        {
                            tempState(currIndex1) = 0.0;
                            tempCovar(currIndex1,currIndex1) = var1.getInitialVariance();
                        }
                    }
                    // old variable
                    else
                    {
                        tempState(currIndex1) = prevState(prevIndex1);
                        tempCovar(currIndex1,currIndex1) = prevCovar(prevIndex1,prevIndex1);
                    }

                    for(int j=i+1; j<numUnk; ++j)
                    {
                        Variable var2( currVarVec[j] );

                        int prevIndex2( var2.getPreviousIndex() );
                        int currIndex2( var2.getCurrentIndex() );

                        if(prevIndex1 == -1 || prevIndex2 == -1)
                        {
                            tempCovar(currIndex1,currIndex2) = 0.0;
                            tempCovar(currIndex2,currIndex1) = tempCovar(currIndex1,currIndex2);
                        }
                        else
                        {
                            tempCovar(currIndex1,currIndex2) = prevCovar(prevIndex1,prevIndex2);
                            tempCovar(currIndex2,currIndex1) = tempCovar(currIndex1,currIndex2);
                        }
                    }
                }

                //////// end of filter initialization ////////


                //////// start of measurement update ////////

                numObs = numObs + 1;
    //            numObs = numObs + 1 + numSat*3;

                // H
                Eigen::MatrixXd H = Eigen::MatrixXd::Zero(numObs,numUnk);

                // omc
                Eigen::VectorXd omc(numObs);
                for(int i=0; i<numObs; ++i) omc(i) = 0.0;

                // noise
                Eigen::MatrixXd R = Eigen::MatrixXd::Zero(numObs,numObs);

                int row(0);
                double obs(0.0), com(0.0);

                obs = 0.0; com = 0.0;
                for(SatIDSet::iterator it = satSet.begin();
                    it != satSet.end();
                    ++it)
                {
                    SatID sat(*it);

                    if(satIndexMap.find(sat) == satIndexMap.end()) continue;

                    int index( satIndexMap[sat] );

                    if(finalSatClock.find(sat) == finalSatClock.end()) continue;

                    obs += finalSatClock[sat];
                    com += tempState(index+6+numSRPC);

                    H(row,index+6+numSRPC) = 1.0;
                }

                omc(row) = obs - com;
                R(row,row) = 1e-1;

                row += 1;

    /*
                obs = 0.0; com = 0.0;
                for(SatIDSet::iterator it = satSet.begin();
                    it != satSet.end();
                    ++it)
                {
                    SatID sat(*it);

                    if(satIndexMap.find(sat) == satIndexMap.end()) continue;

                    int index( satIndexMap[sat] );

                    if(finalSatOrbitECI.find(sat) == finalSatOrbitECI.end()) continue;

                    for(int i=0; i<3; ++i)
                    {
                        obs = finalSatOrbitECI[sat](i);
                        com = tempState(index+i);

                        omc(row+i) = obs - com;
                        H(row+i, index+i) = 1.0;

                        R(row+i,row+i) = 5e-1*5e-1;
                    }

                    row += 3;
                }
    */
                double a12( GAMMA_GPS_L1L2/(GAMMA_GPS_L1L2 - 1.0) );
                double b12( -1.0          /(GAMMA_GPS_L1L2 - 1.0) );
                double X( a12*a12 + b12*b12 );

                double cdtr(0.0), zwdr(0.0);
                double cdts(0.0);
                double ambi(0.0);

                double prefitC(0.0), prefitL(0.0);
                double weightC(0.0), weightL(0.0);
                double postfitC(0.0), postfitL(0.0);

                Matrix<double> I3 = ident<double>(3);

                for(gnssDataMap::iterator gdmIt = gData.begin();
                    gdmIt != gData.end();
                    ++gdmIt)
                {
                    for(sourceDataMap::iterator sdmIt = gdmIt->second.begin();
                        sdmIt != gdmIt->second.end();
                        ++sdmIt)
                    {
                        SourceID source( sdmIt->first );

                        if(sourceIndexMap.find(source) == sourceIndexMap.end()) continue;
                        int rxIndex( sourceIndexMap[source] );

                        cdtr = tempState(rxIndex+0);
                        zwdr = tempState(rxIndex+1);

                        if(sourcePosECI.find(source) == sourcePosECI.end()) continue;
                        Vector<double> posSourceECI( sourcePosECI[source] );

                        if(sourceClock.find(source) == sourceClock.end()) continue;
                        double clock( sourceClock[source] );

                        for(satTypeValueMap::iterator stvmIt = sdmIt->second.begin();
                            stvmIt != sdmIt->second.end();
                            ++stvmIt)
                        {
                            SatID sat( stvmIt->first );
                            typeValueMap tvm( stvmIt->second );

                            if(satIndexMap.find(sat) == satIndexMap.end()) continue;

                            if(tvm.find(TypeID::wetMap) == tvm.end()) continue;
                            double wmf( tvm[TypeID::wetMap] );

                            if(tvm.find(TypeID::rho) == tvm.end()) continue;
                            double dtRho( tvm[TypeID::rho]/C_MPS );

                            double dT( clock/C_MPS + dtRho );
                            double dT2( dT * dT );

                            if( tvm.find(TypeID::prefitC12ForPOD) == tvm.end() ||
                                tvm.find(TypeID::prefitL12ForPOD) == tvm.end() )
                                continue;

                            prefitC = tvm[TypeID::prefitC12ForPOD];
                            prefitL = tvm[TypeID::prefitL12ForPOD];

                            if( tvm.find(TypeID::weight) == tvm.end() )
                                continue;

                            weightC = tvm[TypeID::weight]*1e+0;
                            weightL = tvm[TypeID::weight]*1e+4;

                            int svIndex( satIndexMap[sat] );
                            int rxsvIndex( sourceSatIndexMap[source][sat] );

                            cdts = tempState(svIndex+6+numSRPC);
                            ambi = tempState(rxsvIndex+0);


                            Vector<double> rs(3,0.0), vs(3,0.0);
                            for(int i=0; i<3; ++i)
                            {
                                rs(i) = tempState(svIndex+i+0);
                                vs(i) = tempState(svIndex+i+3);
                            }

                            double r1 = norm(rs);
                            double r2 = r1 * r1;
                            double r3 = r2 * r1;

                            double k = GM_EARTH/r3;

                            Vector<double> posSatECI = rs - vs*dT - 0.5*k*rs*dT2;

                            double rho = norm(posSatECI-posSourceECI);

                            Matrix<double> rrt = outer(rs,rs);

                            Matrix<double> phi = I3 - 0.5*k*(I3-3.0*rrt/r2)*dT2;

                            Vector<double> par = (posSatECI-posSourceECI)/rho;

                            Vector<double> rho_rsat = phi*par;
                            Vector<double> rho_vsat = -par*dT;

                            H(row+0,rxIndex+0) = +1.0;
                            H(row+0,rxIndex+1) =  wmf;
                            for(int i=0; i<3; ++i)
                            {
                                H(row+0,svIndex+i+0) = rho_rsat(i);
                                H(row+0,svIndex+i+3) = rho_vsat(i);
                            }
                            H(row+0,svIndex+6+numSRPC) = -1.0;

                            H(row+1,rxIndex+0) = +1.0;
                            H(row+1,rxIndex+1) =  wmf;
                            for(int i=0; i<3; ++i)
                            {
                                H(row+1,svIndex+i+0) = rho_rsat(i);
                                H(row+1,svIndex+i+3) = rho_vsat(i);
                            }
                            H(row+1,svIndex+6+numSRPC) = -1.0;
                            H(row+1,rxsvIndex+0) = -1.0;

                            omc(row+0) = prefitC - rho - clock - cdtr - wmf*zwdr + cdts;
                            omc(row+1) = prefitL - rho - clock - cdtr - wmf*zwdr + cdts + ambi;

                            R(row+0,row+0) = X/weightC*9e-2;
                            R(row+1,row+1) = X/weightL*9e-2;

                            row += 2;

                        } // End of for(satTypeValueMap::iterator stvmIt = ...)

                    } // End of for(sourceDataMap::iterator sdmIt = ...)

                } // End of for(gnssDataMap::iterator gdmIt = ...)


                Eigen::MatrixXd D( tempCovar*H.transpose() );
                Eigen::MatrixXd S( H*D + R );
                Eigen::MatrixXd U( S.llt().matrixU() );
                Eigen::MatrixXd L( U.transpose() );
                Eigen::MatrixXd E( L.triangularView<Eigen::Lower>().solve(D.transpose()).transpose() );
                Eigen::MatrixXd K( U.triangularView<Eigen::Upper>().solve(E.transpose()).transpose() );

                tempState = tempState + K*omc;
                tempCovar = tempCovar - E*E.transpose();


                // postfit
                gnssDataMap tData( gData );
                bool isValid( true );
                for(gnssDataMap::iterator gdmIt = tData.begin();
                    gdmIt != tData.end();
                    ++gdmIt)
                {
                    for(sourceDataMap::iterator sdmIt = gdmIt->second.begin();
                        sdmIt != gdmIt->second.end();
                        ++sdmIt)
                    {
                        SourceID source( sdmIt->first );

                        if(sourceIndexMap.find(source) == sourceIndexMap.end()) continue;
                        int rxIndex( sourceIndexMap[source] );

                        cdtr = tempState(rxIndex+0);
                        zwdr = tempState(rxIndex+1);

                        if(sourcePosECI.find(source) == sourcePosECI.end()) continue;
                        Vector<double> posSourceECI( sourcePosECI[source] );

                        if(sourceClock.find(source) == sourceClock.end()) continue;
                        double clock( sourceClock[source] );

                        for(satTypeValueMap::iterator stvmIt = sdmIt->second.begin();
                            stvmIt != sdmIt->second.end();
                            ++stvmIt)
                        {
                            SatID sat( stvmIt->first );
                            typeValueMap tvm( stvmIt->second );

                            if(satIndexMap.find(sat) == satIndexMap.end()) continue;

                            if(tvm.find(TypeID::wetMap) == tvm.end()) continue;
                            double wmf( tvm[TypeID::wetMap] );

                            if(tvm.find(TypeID::rho) == tvm.end()) continue;
                            double dtRho( tvm[TypeID::rho]/C_MPS );

                            double dT( clock/C_MPS + dtRho );
                            double dT2( dT * dT );

                            if( tvm.find(TypeID::prefitC12ForPOD) == tvm.end() ||
                                tvm.find(TypeID::prefitL12ForPOD) == tvm.end() )
                                continue;

                            prefitC = tvm[TypeID::prefitC12ForPOD];
                            prefitL = tvm[TypeID::prefitL12ForPOD];

                            if( tvm.find(TypeID::weight) == tvm.end() )
                                continue;

                            weightC = tvm[TypeID::weight]*1e+0;
                            weightL = tvm[TypeID::weight]*1e+4;

                            int svIndex( satIndexMap[sat] );
                            int rxsvIndex( sourceSatIndexMap[source][sat] );

                            cdts = tempState(svIndex+6+numSRPC);
                            ambi = tempState(rxsvIndex+0);


                            Vector<double> rs(3,0.0), vs(3,0.0);
                            for(int i=0; i<3; ++i)
                            {
                                rs(i) = tempState(svIndex+i+0);
                                vs(i) = tempState(svIndex+i+3);
                            }

                            double r1 = norm(rs);
                            double r2 = r1 * r1;
                            double r3 = r2 * r1;

                            double k = GM_EARTH/r3;

                            Vector<double> posSatECI = rs - vs*dT - 0.5*k*rs*dT2;

                            double rho = norm(posSatECI-posSourceECI);

                            postfitC = prefitC - rho - clock - cdtr - wmf*zwdr + cdts;
                            postfitL = prefitL - rho - clock - cdtr - wmf*zwdr + cdts + ambi;

                            if( std::abs(postfitC) > 20.0 ||
                                std::abs(postfitL) > 0.05 )
                            {
                                gData.removeSatID(source,sat);
                                isValid = false;
                                continue;
                            }

                        } // End of for(satTypeValueMap::iterator stvmIt = ...)

                    } // End of for(sourceDataMap::iterator sdmIt = ...)

                } // End of for(gnssDataMap::iterator gdmIt = ...)


                SourceIDSet sourceRejectedSet;

                for(gnssDataMap::iterator gdmIt = gData.begin();
                    gdmIt != gData.end();
                    ++gdmIt)
                {
                    for(sourceDataMap::iterator sdmIt = gdmIt->second.begin();
                        sdmIt != gdmIt->second.end();
                        ++sdmIt)
                    {
                        SourceID source( sdmIt->first );
                        satTypeValueMap stvm( sdmIt->second );

                        if(stvm.size() <= 4) sourceRejectedSet.insert(source);
                    }
                }

                gData.removeSourceID( sourceRejectedSet );

                if( isValid ) break;

            }


            currState = tempState;
            currCovar = tempCovar;


            double clock3( Counter::now() );

            debugStream << setw(4) << CivilTime(gps).year
                        << setw(3) << CivilTime(gps).month
                        << setw(3) << CivilTime(gps).day
                        << setw(3) << CivilTime(gps).hour
                        << setw(3) << CivilTime(gps).minute
                        << setprecision(2)
                        << setw(6) << CivilTime(gps).second;

            debugStream << setw(4) << numSource
                        << setw(4) << numSat
                        << setw(6) << numUnk
                        << setw(6) << numObs;

            debugStream << setw(6) << clock2-clock1
                        << setw(6) << clock3-clock2;

            Vector<double> rcom_c(3,0.0), rcom_t(3,0.0);
            Vector<double> rref_t(3,0.0), rref_c(3,0.0);
            Vector<double> vref_t(3,0.0), vref_c(3,0.0);
            double clk_com(0.0), clk_ref(0.0);

            for(SatIDSet::iterator it = allSatSet.begin();
                it != allSatSet.end();
                ++it)
            {
                SatID sat( *it );

                Vector<double> dXYZ(3,0.0), dRTN(3,0.0);
                Vector<double> srpc(numSRPC,0.0);
                int num(0);

                YawData yawData;

                if(satIndexMap.find(sat) != satIndexMap.end())
                {
                    try
                    {
                        int index(satIndexMap[sat]);

                        num = satSourceNumMap[sat];

                        rcom_c(0) = currState(index+0);
                        rcom_c(1) = currState(index+1);
                        rcom_c(2) = currState(index+2);

                        rref_t = sp3Store.getXvt(sat,gps).x.toVector();
                        vref_t = sp3Store.getXvt(sat,gps).v.toVector();

                        rref_c = t2cRaw * rref_t;
                        vref_c = t2cRaw * vref_t + t2cDot * rref_t;

                        dXYZ = rcom_c - rref_c;

                        dRTN = refSys.XYZ2RTN(dXYZ, rref_c, vref_c);

                        for(int i=0; i<numSRPC; ++i)
                        {
                            srpc(i) = currState(index+i+6);
                        }

                        clk_com = currState(index+6+numSRPC)/C_MPS;
                        clk_ref = sp3Store.getXvt(sat,gps).clkbias;

                        yawData = satYawData[sat];

                        oldSatSRPC[sat] = srpc;
                        oldSatClock[sat] = clk_com;
                    }
                    catch(...)
                    {
                        dXYZ.resize(3,0.0);
                        dRTN.resize(3,0.0);
                        srpc.resize(numSRPC,0.0);
                        clk_com = 0.0;
                        clk_ref = 0.0;
                    }
                }
                else
                {
                    dXYZ.resize(3,0.0);
                    dRTN.resize(3,0.0);
                    srpc.resize(numSRPC,0.0);
                    clk_com = 0.0;
                    clk_ref = 0.0;
                }


                debugStream << setprecision(3);

                debugStream << setw( 4) << sat.id;

                debugStream << setw( 4) << num;

                debugStream << setw(10) << yawData.beta
                            << setw( 4) << yawData.event
                            << setw(10) << yawData.nominal
                            << setw(10) << yawData.modeled;

                debugStream << setw(8) << dRTN(0)
                            << setw(8) << dRTN(1)
                            << setw(8) << dRTN(2);

//                debugStream << setprecision(5);
//                for(int i=0; i<numSRPC; ++i)
//                {
//                    debugStream << setw(12) << srpc(i)*1e+2;
//                }

                debugStream << setw(15) << clk_com*1e+9
                            << setw(15) << clk_ref*1e+9;
            }

            debugStream << endl;

//            break;


            for(SourceIDSet::iterator it = allSourceSet.begin();
                it != allSourceSet.end();
                ++it)
            {
                SourceID source(*it);

                if(sourceSet.find(source) != sourceSet.end())
                {
                    int index = sourceIndexMap[source];

                    double zwdState = currState(index+1);
                    double zwdCovar = currCovar(index+1,index+1);

                    sourceZWDData[source].gps = gps;
                    sourceZWDData[source].state = zwdState;
                    sourceZWDData[source].covar = zwdCovar;
                }
            }

            //////// end of measment update ////////


            //////// start of time update ////////

            satOrbit.clear();
            satSRPC.clear();
            satClock.clear();

            for(SatIDSet::const_iterator it = allSatSet.begin();
                it != allSatSet.end();
                ++it)
            {
                SatID sat(*it);

                Vector<double> orbit(42+6*numSRPC,0.0);
                orbit( 6) = 1.0; orbit(10) = 1.0; orbit(14) = 1.0;
                orbit(33) = 1.0; orbit(37) = 1.0; orbit(41) = 1.0;

                Vector<double> srpc(numSRPC,0.0);
                srpc(0) = -1.0;

                double clock(0.0);

                Vector<double> rsat_t(3,0.0), vsat_t(3,0.0);
                Vector<double> rsat_c(3,0.0), vsat_c(3,0.0);

                if(satIndexMap.find(sat) != satIndexMap.end())
                {
                    int index( satIndexMap[sat] );

                    for(int i=0; i<6; ++i)
                    {
                        orbit(i) = currState(index+i+0);
                    }
                    for(int i=0; i<numSRPC; ++i)
                    {
                        srpc(i) = currState(index+i+6);
                    }

                    clock = currState(index+6+numSRPC);

                    satOrbit[sat] = orbit;
                    satSRPC[sat] = srpc;
                    satClock[sat] = clock;
                }
            }

            srp.setSRPCoeff( satSRPC );
            gnss.setSRPModel( srp );
            rkf78.setEquationOfMotion( gnss );

            rkf78.setCurrentTime( tt );
            rkf78.setCurrentState( satOrbit );

            satOrbit = rkf78.integrateTo( tt+dt );


            // source-related variable
            for(SourceIDSet::iterator it = sourceSet.begin();
                it != sourceSet.end();
                ++it)
            {
                SourceID source( *it );

                int index(sourceIndexMap[source]);

                currState(index+0) = 0.0;

                currCovar(index+0,index+0) = 1e+2*1e+2;
                currCovar(index+1,index+1) += 3e-8*dt;
            }


            // sat-related variable
            for(SatIDSet::iterator it = satSet.begin();
                it != satSet.end();
                ++it)
            {
                SatID sat( *it );

                if(satIndexMap.find(sat) == satIndexMap.end()) continue;

                int index(satIndexMap[sat]);

                Vector<double> orbit( satOrbit[sat] );

                // state
                for(int i=0; i<6; ++i)
                {
                    currState(index+i) = orbit(i);
                }

//                currState(index+6+numSRPC) = 0.0;

                // phi
                int num(6+numSRPC);

                Eigen::MatrixXd phi = Eigen::MatrixXd::Identity(num,num);

                for(int i=0; i<3; ++i)
                {
                    for(int j=0; j<3; ++j)
                    {
                        phi(i+0,j+0) = orbit( 6+3*i+j); // dr1/dr0
                        phi(i+0,j+3) = orbit(15+3*i+j); // dr1/dv0
                        phi(i+3,j+0) = orbit(24+3*i+j); // dv1/dr0
                        phi(i+3,j+3) = orbit(33+3*i+j); // dv1/dv0
                    }

                    for(int j=0; j<numSRPC; ++j)
                    {
                        phi(i+0,j+6) = orbit(42+numSRPC*0+numSRPC*i+j); // dr1/dp0
                        phi(i+3,j+6) = orbit(42+numSRPC*3+numSRPC*i+j); // dv1/dp0
                    }
                }

                // phi * p
                Eigen::MatrixXd pSat1 = Eigen::MatrixXd::Zero(num,numUnk);

                for(int i=0; i<num; ++i)
                {
                    for(int j=0; j<numUnk; ++j)
                    {
                        pSat1(i,j) = currCovar(index+i,j);
                    }
                }

                pSat1 = phi * pSat1;

                for(int i=0; i<num; ++i)
                {
                    for(int j=0; j<numUnk; ++j)
                    {
                        currCovar(index+i,j) = pSat1(i,j);
                    }
                }

                // phi * p * phi'
                Eigen::MatrixXd pSat2 = Eigen::MatrixXd::Zero(numUnk,num);

                for(int i=0; i<numUnk; ++i)
                {
                    for(int j=0; j<num; ++j)
                    {
                        pSat2(i,j) = currCovar(i,index+j);
                    }
                }

                pSat2 = pSat2 * phi.transpose();

                for(int i=0; i<numUnk; ++i)
                {
                    for(int j=0; j<num; ++j)
                    {
                        currCovar(i,index+j) = pSat2(i,j);
                    }
                }

                // phi * p * phi' + q
                YawData yawData( satYawData[sat] );

                int factor(1e+0);
                if(yawData.event != 0) factor = 1e+2;

                for(int i=0; i<numSRPC; ++i)
                {
                    currCovar(index+i+6,index+i+6) += srpcNoise(i)*factor;
                }

                currCovar(index+6+numSRPC,index+6+numSRPC) += 1e+2*1e+2;
            }

            for(int i=0; i<numUnk; ++i)
            {
                for(int j=i+1; j<numUnk; ++j)
                {
                    currCovar(j,i) = currCovar(i,j);
                }
            }

/*
            for(int i=0; i<numUnk; ++i)
            {
                Variable var( currVarVec[i] );
                int currIndex( var.getCurrentIndex() );

                cout << setw(20) << var.getType()
                     << setw(5) << var.getSourceIndexed()
                     << setw(5) << var.getSatIndexed()
                     << setw(20) << currState(currIndex)
                     << setw(20) << currCovar(currIndex,currIndex)
                     << endl;
            }

            break;
*/
            //////// end of time update ////////

            prevVarVec = currVarVec;
            prevState = currState;
            prevCovar = currCovar;

            double clock4( Counter::now() );

            cerr << fixed;

            cerr << setw(4) << CivilTime(gps).year
                 << setw(3) << CivilTime(gps).month
                 << setw(3) << CivilTime(gps).day
                 << setw(3) << CivilTime(gps).hour
                 << setw(3) << CivilTime(gps).minute
                 << setprecision(2)
                 << setw(6) << CivilTime(gps).second;

            cerr << setprecision(2)
                 << setw(12) << "prepare: "
                 << setw(6) << clock2-clock1
                 << setw(15) << "meas update: "
                 << setw(6) << clock3-clock2
                 << setw(15) << "time update: "
                 << setw(6) << clock4-clock3
                 << endl;

            gps += 30.0;

            if(gps >= gps2) break;

        }
        catch(...)
        {
            cerr << "processing error." << endl;
            break;
        }

    } // End of 'while( true )'

    debugStream.close();

    double clock_end( Counter::now() );

    cerr << fixed;
    cerr << setw(20) << "time elapsed: "
         << setw(10) << setprecision(3)
         << clock_end - clock_start << endl;

    return 0;

} // main()
