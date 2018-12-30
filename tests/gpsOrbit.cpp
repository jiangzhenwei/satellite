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
    markCSMW.setSatSystem(SatID::systemGPS);
    markCSMW.setObsType(TypeID::MW12);
    markCSMW.setLLIType1(TypeID::LLI1);
    markCSMW.setLLIType2(TypeID::LLI2);
    markCSMW.setResultType1(TypeID::CSL1);
    markCSMW.setResultType2(TypeID::CSL2);
    markCSMW.setDeltaTMax( 61.0 );
    markCSMW.setMaxNumLambdas( 2.0*WL_WAVELENGTH_GPS_L1L2 );
    markCSMW.setUseLLI( true );


    // LICSDetector
    LICSDetector markCSLI;
    markCSLI.setSatSystem(SatID::systemGPS);
    markCSLI.setObsType(TypeID::LI12);
    markCSLI.setLLIType1(TypeID::LLI1);
    markCSLI.setLLIType2(TypeID::LLI2);
    markCSLI.setResultType1(TypeID::CSL1);
    markCSLI.setResultType2(TypeID::CSL2);
    markCSLI.setDeltaTMax( 61.0 );
    markCSLI.setMinThreshold( 0.240 );
    markCSLI.setLIDrift( 0.002 );
    markCSLI.setUseLLI( true );


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
    staClock.setMaxRMSOfGPS(5.0);


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

    Vector<double> prevState, currState;
    Matrix<double> prevCovar, currCovar;


    // RX related variables
    Variable clkSta(TypeID::dcdtSta);
    clkSta.setSourceIndexed( true );
    clkSta.setSatIndexed( false );
    clkSta.setInitialVariance( 3e+5*3e+5 );

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
    vxSat.setInitialVariance( 1e-3*1e-3 );

    Variable vySat(TypeID::dSatVY);
    vySat.setSourceIndexed( false );
    vySat.setSatIndexed( true );
    vySat.setInitialVariance( 1e-3*1e-3 );

    Variable vzSat(TypeID::dSatVZ);
    vzSat.setSourceIndexed( false );
    vzSat.setSatIndexed( true );
    vzSat.setInitialVariance( 1e-3*1e-3 );

    Variable d0Sat(TypeID::dSRPC1);
    d0Sat.setSourceIndexed( false );
    d0Sat.setSatIndexed( true );
    d0Sat.setInitialVariance( 1e+0*1e+0 );

    Variable dcSat(TypeID::dSRPC2);
    dcSat.setSourceIndexed( false );
    dcSat.setSatIndexed( true );
    dcSat.setInitialVariance( 1e+0*1e+0 );

    Variable dsSat(TypeID::dSRPC3);
    dsSat.setSourceIndexed( false );
    dsSat.setSatIndexed( true );
    dsSat.setInitialVariance( 1e+0*1e+0 );

    Variable y0Sat(TypeID::dSRPC4);
    y0Sat.setSourceIndexed( false );
    y0Sat.setSatIndexed( true );
    y0Sat.setInitialVariance( 1e+0*1e+0 );

    Variable ycSat(TypeID::dSRPC5);
    ycSat.setSourceIndexed( false );
    ycSat.setSatIndexed( true );
    ycSat.setInitialVariance( 1e+0*1e+0 );

    Variable ysSat(TypeID::dSRPC6);
    ysSat.setSourceIndexed( false );
    ysSat.setSatIndexed( true );
    ysSat.setInitialVariance( 1e+0*1e+0 );

    Variable b0Sat(TypeID::dSRPC7);
    b0Sat.setSourceIndexed( false );
    b0Sat.setSatIndexed( true );
    b0Sat.setInitialVariance( 1e+0*1e+0 );

    Variable bcSat(TypeID::dSRPC8);
    bcSat.setSourceIndexed( false );
    bcSat.setSatIndexed( true );
    bcSat.setInitialVariance( 1e+0*1e+0 );

    Variable bsSat(TypeID::dSRPC9);
    bsSat.setSourceIndexed( false );
    bsSat.setSatIndexed( true );
    bsSat.setInitialVariance( 1e+0*1e+0 );

    Variable clkSat(TypeID::dcdtSat);
    clkSat.setSourceIndexed( false );
    clkSat.setSatIndexed( true );
    clkSat.setInitialVariance( 3e+5*3e+5 );

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
    double srpcNoise;

    try
    {
        srpcNoise = confReader.getValueAsDouble("srpcNoise", "DEFAULT");
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

    ofstream debugStrm( debugFileName.c_str() );

    debugStrm << fixed;

    cout << fixed;

    map<SourceID, varData> sourceTropoData;

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
        double sow( gws.sow );
        int dow( int(sow/86400.0) );

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
                cerr << "sp3/clk file load error." << endl;
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
                cerr << "dcb file open error." << endl;
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

                    cerr << "obs file read error." << endl;
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

//                cout << setw(5) << source.sourceName;
//                cout << setw(30) << CivilTime(gRin.header.epoch);
//                cout << setw(5) << gRin.body.numSats();
//                cout << endl;

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

//                cout << setw(5) << source.sourceName;
//                cout << setw(30) << CivilTime(gRin.header.epoch);
//                cout << setw(5) << gRin.body.numSats();
//                cout << endl;

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

            currVarVec.clear();

            sourceIndexMap.clear();
            satIndexMap.clear();
            sourceSatIndexMap.clear();
/*
            for(gnssDataMap::iterator gdmIt = gData.begin();
                gdmIt != gData.end();
                ++gdmIt)
            {
                for(sourceDataMap::iterator sdmIt = gdmIt->second.begin();
                    sdmIt != gdmIt->second.end();
                    ++sdmIt)
                {
                    cout << setw(5) << sdmIt->first.sourceName << endl;

                    for(satTypeValueMap::iterator stvmIt = sdmIt->second.begin();
                        stvmIt != sdmIt->second.end();
                        ++stvmIt)
                    {
                        cout << setw(5) << stvmIt->first << endl;

                        for(typeValueMap::iterator tvmIt = stvmIt->second.begin();
                            tvmIt != stvmIt->second.end();
                            ++tvmIt)
                        {
                            cout << setw(15) << tvmIt->first
                                 << setw(20) << tvmIt->second
                                 << endl;

                        }
                    }

                    break;
                }

                break;
            }

            break;
*/

            //////// start of filter initialization ////////

            // number of source
            sourceSet = gData.getSourceIDSet();
            numSource = sourceSet.size();

            // number of sat
            satSet = gData.getSatIDSet();
            numSat = satSet.size();

            if(numSource == 0 || numSat == 0) break;


            map<SatID,int> satSourceNumMap;

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
            int numUnk( currVarVec.size() );
            currState.resize(numUnk,0.0);
            currCovar.resize(numUnk,numUnk,0.0);

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
                            currState(currIndex1) = finalSatOrbitECI[sat](0);
                        }
                        currCovar(currIndex1,currIndex1) = var1.getInitialVariance();
                    }
                    else if(type == TypeID::dSatY)
                    {
                        SatID sat( var1.getSatellite() );

                        if(finalSatOrbitECI.find(sat) != finalSatOrbitECI.end())
                        {
                            currState(currIndex1) = finalSatOrbitECI[sat](1);
                        }
                        currCovar(currIndex1,currIndex1) = var1.getInitialVariance();
                    }
                    else if(type == TypeID::dSatZ)
                    {
                        SatID sat( var1.getSatellite() );

                        if(finalSatOrbitECI.find(sat) != finalSatOrbitECI.end())
                        {
                            currState(currIndex1) = finalSatOrbitECI[sat](2);
                        }
                        currCovar(currIndex1,currIndex1) = var1.getInitialVariance();
                    }
                    else if(type == TypeID::dSatVX)
                    {
                        SatID sat( var1.getSatellite() );

                        if(finalSatOrbitECI.find(sat) != finalSatOrbitECI.end())
                        {
                            currState(currIndex1) = finalSatOrbitECI[sat](3);
                        }
                        currCovar(currIndex1,currIndex1) = var1.getInitialVariance();
                    }
                    else if(type == TypeID::dSatVY)
                    {
                        SatID sat( var1.getSatellite() );

                        if(finalSatOrbitECI.find(sat) != finalSatOrbitECI.end())
                        {
                            currState(currIndex1) = finalSatOrbitECI[sat](4);
                        }
                        currCovar(currIndex1,currIndex1) = var1.getInitialVariance();
                    }
                    else if(type == TypeID::dSatVZ)
                    {
                        SatID sat( var1.getSatellite() );

                        if(finalSatOrbitECI.find(sat) != finalSatOrbitECI.end())
                        {
                            currState(currIndex1) = finalSatOrbitECI[sat](5);
                        }
                        currCovar(currIndex1,currIndex1) = var1.getInitialVariance();
                    }
                    else if(type == TypeID::dcdtSat)
                    {
                        SatID sat( var1.getSatellite() );

                        if(finalSatClock.find(sat) != finalSatClock.end())
                        {
                            currState(currIndex1) = finalSatClock[sat];
                            currCovar(currIndex1,currIndex1) = 1e+2*1e+2;
                        }
                        else
                        {
                            currState(currIndex1) = 0.0;
                            currCovar(currIndex1,currIndex1) = 3e+5*3e+5;
                        }
                    }
                    else if(type == TypeID::dSRPC1)
                    {
                        SatID sat( var1.getSatellite() );

                        currState(currIndex1) = -1.0;
                        currCovar(currIndex1,currIndex1) = var1.getInitialVariance();
                    }
                    else if(type == TypeID::wetMap)
                    {
                        SourceID source( var1.getSource() );

                        if(sourceTropoData.find(source) != sourceTropoData.end())
                        {
                            CommonTime gpsZWD(CommonTime::BEGINNING_OF_TIME);

                            double stateZWD(0.0), covarZWD(0.0);

                            gpsZWD = sourceTropoData[source].gps;

                            stateZWD = sourceTropoData[source].state;
                            covarZWD = sourceTropoData[source].covar;

                            covarZWD += 3e-8*(gps-gpsZWD);

                            if(covarZWD > 5e-1*5e-1)
                            {
                                stateZWD = 0.0;
                                covarZWD = 5e-1*5e-1;
                            }

                            currState(currIndex1) = stateZWD;
                            currCovar(currIndex1,currIndex1) = covarZWD;
                        }
                        else
                        {
                            currState(currIndex1) = 0.0;
                            currCovar(currIndex1,currIndex1) = var1.getInitialVariance();
                        }
                    }
                    else
                    {
                        currState(currIndex1) = 0.0;
                        currCovar(currIndex1,currIndex1) = var1.getInitialVariance();
                    }
                }
                // old variable
                else
                {
                    currState(currIndex1) = prevState(prevIndex1);
                    currCovar(currIndex1,currIndex1) = prevCovar(prevIndex1,prevIndex1);
                }

                for(int j=i+1; j<numUnk; ++j)
                {
                    Variable var2( currVarVec[j] );

                    int prevIndex2( var2.getPreviousIndex() );
                    int currIndex2( var2.getCurrentIndex() );

                    if(prevIndex1 == -1 || prevIndex2 == -1)
                    {
                        currCovar(currIndex1,currIndex2) = 0.0;
                        currCovar(currIndex2,currIndex1) = currCovar(currIndex1,currIndex2);
                    }
                    else
                    {
                        currCovar(currIndex1,currIndex2) = prevCovar(prevIndex1,prevIndex2);
                        currCovar(currIndex2,currIndex1) = currCovar(currIndex1,currIndex2);
                    }
                }
            }

            double clock3( Counter::now() );
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
            //////// end of filter initialization ////////


            //////// start of clock constraint ////////

            double obs(0.0), com(0.0);

            // omc, weight
            double omc(0.0), weight(0.0);

            // h
            Vector<double> h(numUnk,0.0);

            // p * h'
            Vector<double> pht(numUnk,0.0);

            // h * p * h'
            double hpht(0.0);

            // beta
            double beta(0.0);

            // kalman gain
            Vector<double> gamma(numUnk,0.0);

            for(SatIDSet::iterator it = satSet.begin();
                it != satSet.end();
                ++it)
            {
                SatID sat(*it);

                if(satIndexMap.find(sat) == satIndexMap.end()) continue;

                int index( satIndexMap[sat] );

                if(finalSatClock.find(sat) == finalSatClock.end()) continue;

                obs += finalSatClock[sat];
                com += currState(index+6+numSRPC);

                h(index+6+numSRPC) = 1.0;
            }

            omc = obs - com;

            for(int i=0; i<numUnk; ++i)
            {
                for(int j=0; j<numUnk; ++j)
                {
                    if(h(j) != 0.0) pht(i) += currCovar(i,j) * h(j);
                }
            }

            for(int i=0; i<numUnk; ++i)
            {
                if(h(i) != 0.0) hpht += h(i) * pht(i);
            }

            weight = 1e-1;

            beta = 1.0/weight + hpht;

            gamma = pht/beta;

            // state update
            currState = currState + gamma*omc;

            // covariance update
#pragma omp parallel for num_threads(numThreads)
            for(int i=0; i<numUnk; ++i)
            {
                currCovar(i,i) = currCovar(i,i) - gamma(i)*pht(i);

                for(int j=i+1; j<numUnk; ++j)
                {
                    currCovar(i,j) = currCovar(i,j) - gamma(i)*pht(j);
                    currCovar(j,i) = currCovar(i,j);
                }
            }

            //////// end of clock constraint ////////


            //////// start of measurement update ////////

            double cdtr(0.0), zwdr(0.0);
            Vector<double> rSat(3,0.0), vSat(3,0.0);
            double cdts(0.0);
            double ambi(0.0);
            double prefitC(0.0), prefitL(0.0);
            double weightC(0.0), weightL(0.0);

            double r1(0.0), r2(0.0), r3(0.0), k(0.0);
            double rho(0.0);
            Vector<double> rho_rsat(3,0.0), rho_vsat(3,0.0);

            Matrix<double> I = ident<double>(3);
            Matrix<double> rrt(3,3,0.0), phi(3,3,0.0);
            Vector<double> par(3,0.0);

            Vector<double> posSatECI(3,0.0);

            int numObs(0);

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
                        int svIndex( satIndexMap[sat] );

                        int rxsvIndex( sourceSatIndexMap[source][sat] );

                        if(tvm.find(TypeID::wetMap) == tvm.end()) continue;
                        double wmf( tvm[TypeID::wetMap] );

                        double dtSta( clock/C_MPS );

                        if(tvm.find(TypeID::rho) == tvm.end()) continue;
                        double dtRho( tvm[TypeID::rho]/C_MPS );

                        double dT( dtSta + dtRho );
                        double dT2( dT * dT );

                        if(tvm.find(TypeID::prefitL12ForPOD) == tvm.end()) continue;
                        prefitL = tvm[TypeID::prefitL12ForPOD];

                        if(tvm.find(TypeID::prefitC12ForPOD) == tvm.end()) continue;
                        prefitC = tvm[TypeID::prefitC12ForPOD];

                        if(tvm.find(TypeID::weight) == tvm.end()) continue;
                        weightL = tvm[TypeID::weight]*1e+4;
                        weightC = tvm[TypeID::weight]*1e+0;


                        //////// phase ////////

                        cdtr = currState(rxIndex+0);
                        zwdr = currState(rxIndex+1);

                        // h
                        h.resize(numUnk,0.0);
                        h(rxIndex+0) = +1.0;
                        h(rxIndex+1) =  wmf;

                        rSat(0) = currState(svIndex+0);
                        rSat(1) = currState(svIndex+1);
                        rSat(2) = currState(svIndex+2);

                        vSat(0) = currState(svIndex+3);
                        vSat(1) = currState(svIndex+4);
                        vSat(2) = currState(svIndex+5);

                        cdts = currState(svIndex+6+numSRPC);

                        r1 = norm(rSat);
                        r2 = r1 * r1;
                        r3 = r2 * r1;

                        k = GM_EARTH/r3;

                        posSatECI = rSat - vSat*dT - 0.5*k*rSat*dT2;

                        rho = norm(posSatECI-posSourceECI);

                        rrt = outer(rSat,rSat);

                        phi = I - 0.5*k*(I-3.0*rrt/r2)*dT2;

                        par = (posSatECI-posSourceECI)/rho;

                        rho_rsat = phi*par;
                        rho_vsat = -par*dT;

                        h(svIndex+0) = rho_rsat(0);
                        h(svIndex+1) = rho_rsat(1);
                        h(svIndex+2) = rho_rsat(2);
                        h(svIndex+3) = rho_vsat(0);
                        h(svIndex+4) = rho_vsat(1);
                        h(svIndex+5) = rho_vsat(2);
                        h(svIndex+6+numSRPC) = -1.0;

                        ambi = currState(rxsvIndex+0);

                        h(rxsvIndex+0) = -1.0;

                        // omc
                        omc = prefitL - rho - cdtr - wmf*zwdr + cdts + ambi;

                        // p * h'
                        pht.resize(numUnk,0.0);
                        for(int i=0; i<numUnk; ++i)
                        {
                            for(int j=0; j<numUnk; ++j)
                            {
                                if(h(j) != 0.0) pht(i) += currCovar(i,j) * h(j);
                            }
                        }

                        // h * p * h'
                        hpht = 0.0;
                        for(int i=0; i<numUnk; ++i)
                        {
                            if(h(i) != 0.0) hpht += h(i) * pht(i);
                        }

                        beta = 1.0/weightL + hpht;

                        // kalman gain
                        gamma = pht/beta;

                        // state update
                        currState = currState + gamma*omc;

                        // covariance update
#pragma omp parallel for num_threads(numThreads)
                        for(int i=0; i<numUnk; ++i)
                        {
                            currCovar(i,i) = currCovar(i,i) - gamma(i)*pht(i);

                            for(int j=i+1; j<numUnk; ++j)
                            {
                                currCovar(i,j) = currCovar(i,j) - gamma(i)*pht(j);
                                currCovar(j,i) = currCovar(i,j);
                            }
                        }


                        //////// code ////////

                        cdtr = currState(rxIndex+0);
                        zwdr = currState(rxIndex+1);

                        // h
                        h.resize(numUnk,0.0);
                        h(rxIndex+0) = +1.0;
                        h(rxIndex+1) =  wmf;

                        rSat(0) = currState(svIndex+0);
                        rSat(1) = currState(svIndex+1);
                        rSat(2) = currState(svIndex+2);

                        vSat(0) = currState(svIndex+3);
                        vSat(1) = currState(svIndex+4);
                        vSat(2) = currState(svIndex+5);

                        cdts = currState(svIndex+6+numSRPC);

                        r1 = norm(rSat);
                        r2 = r1 * r1;
                        r3 = r2 * r1;

                        k = GM_EARTH/r3;

                        posSatECI = rSat - vSat*dT - 0.5*k*rSat*dT2;

                        rho = norm(posSatECI-posSourceECI);

                        rrt = outer(rSat,rSat);

                        phi = I - 0.5*k*(I-3.0*rrt/r2)*dT2;

                        par = (posSatECI-posSourceECI)/rho;

                        rho_rsat = phi*par;
                        rho_vsat = -par*dT;

                        h(svIndex+0) = rho_rsat(0);
                        h(svIndex+1) = rho_rsat(1);
                        h(svIndex+2) = rho_rsat(2);
                        h(svIndex+3) = rho_vsat(0);
                        h(svIndex+4) = rho_vsat(1);
                        h(svIndex+5) = rho_vsat(2);
                        h(svIndex+6+numSRPC) = -1.0;

                        // omc
                        omc = prefitC - rho - cdtr - wmf*zwdr + cdts;

                        // p * h'
                        pht.resize(numUnk,0.0);
                        for(int i=0; i<numUnk; ++i)
                        {
                            for(int j=0; j<numUnk; ++j)
                            {
                                if(h(j) != 0.0) pht(i) += currCovar(i,j) * h(j);
                            }
                        }

                        // h * p * h'
                        hpht = 0.0;
                        for(int i=0; i<numUnk; ++i)
                        {
                            if(h(i) != 0.0) hpht += h(i) * pht(i);
                        }

                        beta = 1.0/weightC + hpht;

                        // kalman gain
                        gamma = pht/beta;

                        // state update
                        currState = currState + gamma*omc;

                        // covariance update
#pragma omp parallel for num_threads(numThreads)
                        for(int i=0; i<numUnk; ++i)
                        {
                            currCovar(i,i) = currCovar(i,i) - gamma(i)*pht(i);

                            for(int j=i+1; j<numUnk; ++j)
                            {
                                currCovar(i,j) = currCovar(i,j) - gamma(i)*pht(j);
                                currCovar(j,i) = currCovar(i,j);
                            }
                        }

                        numObs += 2;

                    } // End of for(satTypeValueMap::iterator stvmIt = ...)

                } // End of for(sourceDataMap::iterator sdmIt = ...)

            } // End of for(gnssDataMap::iterator gdmIt = ...)


            double clock4( Counter::now() );

            debugStrm   << setw(4) << CivilTime(gps).year
                        << setw(4) << CivilTime(gps).month
                        << setw(4) << CivilTime(gps).day
                        << setw(4) << CivilTime(gps).hour
                        << setw(4) << CivilTime(gps).minute
                        << setprecision(3)
                        << setw(10) << CivilTime(gps).second;

            debugStrm   << setw(5) << numSource
                        << setw(5) << numSat
                        << setw(8) << numUnk
                        << setw(8) << numObs;

            debugStrm   << setw(8) << clock2-clock1
                        << setw(8) << clock4-clock2;

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

                stringstream satstr;
                satstr << "G";
                if( sat.id < 10 ) satstr << "0";
                satstr << sat.id;

                debugStrm   << setprecision(3)
//                            << setw( 4) << satstr.str();
                            << setw( 4) << sat.id;

                debugStrm   << setw( 4) << num;

                debugStrm   << setw(10) << yawData.beta
                            << setw( 4) << yawData.event;
//                            << setw(10) << yawData.nominal
//                            << setw(10) << yawData.modeled;

                debugStrm   << setw(10) << dRTN(0)
                            << setw(10) << dRTN(1)
                            << setw(10) << dRTN(2);

                debugStrm   << setw(15) << clk_com*1e+9
                            << setw(15) << clk_ref*1e+9;
            }

            debugStrm << endl;

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

                    sourceTropoData[source].gps = gps;
                    sourceTropoData[source].state = zwdState;
                    sourceTropoData[source].covar = zwdCovar;
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

            double clock5( Counter::now() );


            // source-related variable
            for(SourceIDSet::iterator it = sourceSet.begin();
                it != sourceSet.end();
                ++it)
            {
                SourceID source( *it );

                int index(sourceIndexMap[source]);

                currState(index+0) = 0.0;

                currCovar(index+0,index+0) = 3e+5*3e+5;
                currCovar(index+1,index+1) += 3e-8*dt;
            }


            // sat-related variable
            for(SatIDSet::iterator it = satSet.begin();
                it != satSet.end();
                ++it)
            {
                SatID sat( *it );

                if(satIndexMap.find(sat) == satIndexMap.end()) continue;

                if(satYawData.find(sat) == satYawData.end()) continue;

                int index(satIndexMap[sat]);

                YawData yawData(satYawData[sat]);

                Vector<double> orbit(42+6*numSRPC,0.0);
                orbit( 6) = 1.0; orbit(10) = 1.0; orbit(14) = 1.0;
                orbit(33) = 1.0; orbit(37) = 1.0; orbit(41) = 1.0;

                orbit = satOrbit[sat];

                // state
                for(int i=0; i<6; ++i)
                {
                    currState(index+i) = orbit(i);
                }

                currState(index+6+numSRPC) = 0.0;


                // phi
                int num(6+numSRPC);

                phi = ident<double>(num);

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
                Matrix<double> pSat1(num,numUnk,0.0);

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
                Matrix<double> pSat2(numUnk,num,0.0);

                for(int i=0; i<numUnk; ++i)
                {
                    for(int j=0; j<num; ++j)
                    {
                        pSat2(i,j) = currCovar(i,index+j);
                    }
                }

                pSat2 = pSat2 * transpose(phi);

                for(int i=0; i<numUnk; ++i)
                {
                    for(int j=0; j<num; ++j)
                    {
                        currCovar(i,index+j) = pSat2(i,j);
                    }
                }

                // phi * p * phi' + q
                for(int i=0; i<numSRPC; ++i)
                {
                    currCovar(index+i+6,index+i+6) += srpcNoise;
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

//            cout << currCovar << endl;

            //////// end of time update ////////

            prevVarVec = currVarVec;
            prevState = currState;
            prevCovar = currCovar;

            double clock6( Counter::now() );

            cerr << fixed << setprecision(3)
                 << setw(15) << "prepare: "
                 << setw(10) << clock3-clock2
                 << setw(15) << "meas update: "
                 << setw(10) << clock4-clock3
                 << setw(15) << "time update: "
                 << setw(10) << clock6-clock4
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

    debugStrm.close();

    double clock_end( Counter::now() );

    cerr << fixed;
    cerr << setw(20) << "time elapsed: "
         << setw(10) << setprecision(3)
         << clock_end - clock_start << endl;

    return 0;

} // main()
