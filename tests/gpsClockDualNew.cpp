#pragma ident "$Id$"

//#include <complex>
//
//#define EIGEN_SUPERLU_SUPPORT
//#define EIGEN_USE_BLAS
//#define EIGEN_USE_LAPACKE
//
//#define LAPACK_COMPLEX_CUSTOM
//#define lapack_complex_float std::complex<float>
//#define lapack_complex_double std::complex<double>

#include "ConfDataReader.hpp"

#include "DataStructures.hpp"

#include "EOPDataStore2.hpp"

#include "LeapSecStore.hpp"

#include "ReferenceSystem.hpp"

#include "SP3EphemerisStore.hpp"

#include "Rinex3EphemerisStore2.hpp"

#include "SolarSystem.hpp"

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

#include "MWFilter.hpp"

#include "BasicModelForPCE.hpp"

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

#include "Epoch.hpp"

#include "Counter.hpp"

#include "Variable.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>


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
    string confFileName("gps_clock.conf");

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


    // station list
    SourceIDSet allSourceSet;

    SourceIDSet atomicClockSet;

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

        string line, station, atomic;

        while( !stationList.eof() && stationList.good() )
        {
            getline(stationList, line);

            if( stationList.eof() ) break;
            if( stationList.bad() ) break;

            if(line.size() != 0)
            {
                station = line.substr(0,4);
                atomic = line.substr(line.length()-1,1);

                SourceID source;

                transform(station.begin(),station.end(),station.begin(),::toupper);

                source = SourceID(SourceID::Mixed, station);

                if(atomic == "Y") atomicClockSet.insert(source);

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
    SatIDSet allSatSet;

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


    // SP3EphemerisStore
    SP3EphemerisStore sp3Store;


    // Rinex3EphemerisStore2
    Rinex3EphemerisStore2 bceStore;


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
    linearCS.setCheckExist(true);
    linearCS.addLinear(SatID::systemGPS, linearComb.mw12CombOfGPS);
    linearCS.addLinear(SatID::systemGPS, linearComb.li12CombOfGPS);

    // ComputeLinear, Phase Code Alignment
    ComputeLinear linearAlign;
    linearAlign.setCheckExist(true);
    linearAlign.addLinear(SatID::systemGPS, linearComb.q1CombOfGPSL1L2);
    linearAlign.addLinear(SatID::systemGPS, linearComb.q2CombOfGPSL1L2);


    // ComputeLinear, PC
    ComputeLinear linearPC;
    linearPC.setCheckExist(true);
    linearPC.addLinear(SatID::systemGPS, linearComb.pc12CombOfGPS);

    // ComputeLinear, LC
    ComputeLinear linearLC;
    linearLC.setCheckExist(true);
    linearLC.addLinear(SatID::systemGPS, linearComb.lc12CombOfGPS);


    // ComputeLinear, prefit PC
    ComputeLinear prefitPC;
    prefitPC.setCheckExist(true);
    prefitPC.addLinear(SatID::systemGPS, linearComb.pc12PrefitOfGPS);


    // ComputeLinear, prefit PC/LC for PCE
    ComputeLinear prefitForPCE;
    prefitForPCE.setCheckExist(true);
    prefitForPCE.addLinear(SatID::systemGPS, linearComb.pc12PrefitOfGPSForPCE);
    prefitForPCE.addLinear(SatID::systemGPS, linearComb.lc12PrefitOfGPSForPCE);


    // MWCSDetector
    MWCSDetector markCSMW;
    markCSMW.setCheckExist(true);
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
    markCSLI.setCheckExist(true);
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


    // MWFilter
    MWFilter filterMW;


    // BasicModelForPCE
    BasicModelForPCE basicModel;
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
    PhaseCodeAlignment phaseAlignGPSL1;
    phaseAlignGPSL1.setSatSystem(SatID::systemGPS);
    phaseAlignGPSL1.setCodeType(TypeID::Q1OfL1L2);
    phaseAlignGPSL1.setPhaseType(TypeID::L1);
    phaseAlignGPSL1.setPhaseWavelength(L1_WAVELENGTH_GPS);
    phaseAlignGPSL1.setUseSatArc(true);

    PhaseCodeAlignment phaseAlignGPSL2;
    phaseAlignGPSL2.setSatSystem(SatID::systemGPS);
    phaseAlignGPSL2.setCodeType(TypeID::Q2OfL1L2);
    phaseAlignGPSL2.setPhaseType(TypeID::L2);
    phaseAlignGPSL2.setPhaseWavelength(L2_WAVELENGTH_GPS);
    phaseAlignGPSL2.setUseSatArc(true);


    SourceIDSet sourceSet;
    SatIDSet satSet;

    int numSource(0), numSat(0);


    VariableSet prevVarSet, currVarSet;

    Eigen::VectorXd prevState, currState;
    Eigen::MatrixXd prevCovar, currCovar;


    VariableSet wlVarSet, nlVarSet;

    Eigen::VectorXd wlState, nlState;
    Eigen::MatrixXd wlCovar, nlCovar;


    // RX related variables
    Variable clkStaVar(TypeID::dcdtSta);
    clkStaVar.setSourceIndexed( true );
    clkStaVar.setSatIndexed( false );
    clkStaVar.setInitialVariance( 1e+2*1e+2 );

    Variable zwdStaVar(TypeID::wetMap);
    zwdStaVar.setSourceIndexed( true );
    zwdStaVar.setSatIndexed( false );
    zwdStaVar.setInitialVariance( 5e-1*5e-1 );


    // SV related variables
    Variable clkSatVar(TypeID::dcdtSat);
    clkSatVar.setSourceIndexed( false );
    clkSatVar.setSatIndexed( true );
    clkSatVar.setInitialVariance( 3e+5*3e+5 );


    // RX-SV related variables
    Variable ambi12Var(TypeID::BLC12);
    ambi12Var.setSourceIndexed( true );
    ambi12Var.setSatIndexed( true );
    ambi12Var.setInitialVariance( 1e+2*1e+2 );


    // UPD
    Variable updStaVar(TypeID::updStaWL12);
    updStaVar.setSourceIndexed( true );
    updStaVar.setSatIndexed( false );

    Variable updSatVar(TypeID::updSatWL12);
    updSatVar.setSourceIndexed( false );
    updSatVar.setSatIndexed( true );


    TypeIDSet keepTypes;
    keepTypes.insert(TypeID::prefitC12ForPCE);
    keepTypes.insert(TypeID::prefitL12ForPCE);
    keepTypes.insert(TypeID::weight);
    keepTypes.insert(TypeID::wetMap);
    keepTypes.insert(TypeID::CSL1);
    keepTypes.insert(TypeID::meanMW12);
    keepTypes.insert(TypeID::varMW12);


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

    // begin time
    year = asInt( time1.substr( 0,4) );
    mon  = asInt( time1.substr( 5,2) );
    day  = asInt( time1.substr( 8,2) );
    hour = asInt( time1.substr(11,2) );
    min  = asInt( time1.substr(14,2) );
    sec  = asDouble( time1.substr(17,2) );

    CivilTime civil1(year,mon,day,hour,min,sec, TimeSystem::GPS);
    CommonTime gps1( civil1.convertToCommonTime() );

    cerr << "begin time: " << civil1 << endl;

    // break time
    year = asInt( time2.substr( 0,4) );
    mon  = asInt( time2.substr( 5,2) );
    day  = asInt( time2.substr( 8,2) );
    hour = asInt( time2.substr(11,2) );
    min  = asInt( time2.substr(14,2) );
    sec  = asDouble( time2.substr(17,2) );

    CivilTime civil2(year,mon,day,hour,min,sec, TimeSystem::GPS);
    CommonTime gps2( civil2.convertToCommonTime() );

    cerr << "break time: " << civil2 << endl;


    // sp3 file name
    string sp3FileName;

    try
    {
        sp3FileName = confReader.getValue("sp3FileName", "DEFAULT");
    }
    catch(...)
    {
        cerr << "sp3 file name get error." << endl;
        exit(-1);
    }
/*
    try
    {
        sp3Store.loadSP3File( sp3FileName );

        basicModel.setSP3Store( sp3Store );
        satPCenter.setEphStore( sp3Store );
        windUp.setEphStore( sp3Store );
    }
    catch(...)
    {
        cerr << "sp3 file load error." << endl;
        exit(-1);
    }
*/

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

    // sp3 file suffix
    string sp3FileSuf;

    try
    {
        sp3FileSuf = confReader.getValue("sp3FileSuf", "DEFAULT");
    }
    catch(...)
    {
        cerr << "sp3 file suffix get error." << endl;
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


    CommonTime gps( gps1 ), utc, tt;

    ofstream outFile("clock_gps_new.txt");

    outFile << fixed;

    cout << fixed;

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

            string sp3FileHdr(sp3FileDir + productAbbr);

            string sp3File, clkFile;

            try
            {
                CommonTime temp;

                // current day - 1
                temp = gps;
                temp -= 86400.0;

                gws = GPSWeekSecond(temp);

                sp3File = sp3FileHdr + asString(gws.week)
                        + asString(gws.getDayOfWeek()) + "." + sp3FileSuf;

                sp3Store.loadSP3File( sp3File );

                // current day
                temp = gps;

                gws = GPSWeekSecond(temp);

                sp3File = sp3FileHdr + asString(gws.week)
                        + asString(gws.getDayOfWeek()) + "." + sp3FileSuf;

                sp3Store.loadSP3File( sp3File );

                // current day + 1
                temp = gps;
                temp += 86400.0;

                gws = GPSWeekSecond(temp);

                sp3File = sp3FileHdr + asString(gws.week)
                        + asString(gws.getDayOfWeek()) + "." + sp3FileSuf;

                sp3Store.loadSP3File( sp3File );
            }
            catch(...)
            {
                cerr << "sp3 and/or clk file load error." << endl;
                break;
            }

            basicModel.setSP3Store( sp3Store );
            satPCenter.setEphStore( sp3Store );
            windUp.setEphStore( sp3Store );

            updateSP3 = false;
        }


        // update NAV, every day
        if(hour==0 && minute==0 && second==0.0) updateNAV = true;

        if( updateNAV )
        {
            bceStore = Rinex3EphemerisStore2();

            string navFileHdr(navFileDir + "brdm");

            string navFile;

            try
            {
                CommonTime temp;
                string dy3;

                // current day - 1
                temp = gps;
                temp -= 86400.0;

                yds = YDSTime(temp);

                if(yds.doy < 10)        dy3 = "00" + asString(yds.doy);
                else if(yds.doy < 100)  dy3 = "0" + asString(yds.doy);
                else                    dy3 = asString(yds.doy);

                navFile = navFileHdr + dy3 + "0." + asString(yds.year).substr(2,2) + "p";

                bceStore.loadFile( navFile );

                // current day
                temp = gps;

                yds = YDSTime(temp);

                if(yds.doy < 10)        dy3 = "00" + asString(yds.doy);
                else if(yds.doy < 100)  dy3 = "0" + asString(yds.doy);
                else                    dy3 = asString(yds.doy);

                navFile = navFileHdr + dy3 + "0." + asString(yds.year).substr(2,2) + "p";

                bceStore.loadFile( navFile );

                // current day + 1
                temp = gps;
                temp += 86400.0;

                yds = YDSTime(temp);

                if(yds.doy < 10)        dy3 = "00" + asString(yds.doy);
                else if(yds.doy < 100)  dy3 = "0" + asString(yds.doy);
                else                    dy3 = asString(yds.doy);

                navFile = navFileHdr + dy3 + "0." + asString(yds.year).substr(2,2) + "p";

                bceStore.loadFile( navFile );
            }
            catch(...)
            {
                cerr << "nav file load error." << endl;
                break;
            }

            basicModel.setBCEStore( bceStore );

            updateNAV = false;
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
                cerr << "msc file load error!" << endl;
                break;
            }

            basicModel.setMSCStore( mscStore );
            computeTM.setMSCStore( mscStore );
            correctObs.setMSCStore( mscStore );
            gravDelay.setMSCStore( mscStore );
            satPCenter.setMSCStore( mscStore );
            windUp.setMSCStore( mscStore );

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

                string yr4( asString(yds.year) );
                string yr2( asString(yds.year).substr(2,2) );
                string dy3;

                if(yds.doy < 10)        dy3 = "00" + asString(yds.doy);
                else if(yds.doy < 100)  dy3 = "0" + asString(yds.doy);
                else                    dy3 = asString(yds.doy);

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
                        (*pObsStream).close();
                        cerr << "obs file '" << file << "' header read error." << endl;
                        continue;
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
                    (*pObsStream).close();
                    cerr << "obs file '" << file << "' data read error." << endl;
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


        gData.keepOnlySatID( allSatSet );

        if(gData.size() == 0) break;


        try
        {
            double clock1( Counter::now() );

            utc = refSys.GPS2UTC( gps );

            Matrix<double> t2cRaw( refSys.T2CMatrix(utc) );
            Matrix<double> t2cDot( refSys.dT2CMatrix(utc) );

            SatIDSet satRejectedSet;

            // sat orbit
            satVectorMap satOrbitECF, satOrbitECI;

            for(SatIDSet::iterator it = allSatSet.begin();
                it != allSatSet.end();
                ++it)
            {
                SatID sat( *it );

                Vector<double> rsat_t(3,0.0), vsat_t(3,0.0);
                Vector<double> rsat_c(3,0.0), vsat_c(3,0.0);
                Vector<double> orbit(6,0.0);

                try
                {
                    rsat_t = sp3Store.getXvt(sat,gps).x.toVector();
                    vsat_t = sp3Store.getXvt(sat,gps).v.toVector();

                    for(int i=0; i<3; ++i)
                    {
                        orbit(i+0) = rsat_t(i);
                        orbit(i+3) = vsat_t(i);
                    }

                    satOrbitECF[sat] = orbit;

                    rsat_c = t2cRaw * rsat_t;
                    vsat_c = t2cRaw * vsat_t + t2cDot * rsat_t;

                    for(int i=0; i<3; ++i)
                    {
                        orbit(i+0) = rsat_c(i);
                        orbit(i+3) = vsat_c(i);
                    }

                    satOrbitECI[sat] = orbit;
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
                  >> prefitPC
                  >> staClock;

            sourceValueMap sourceClock;
            sourceClock = staClock.getSourceClockMap();

            gData >> linearCS               // Cycle Slip
                  >> markCSMW               // MW detector
                  >> markCSLI               // LI detector
                  >> markArc                // obs arc
                  >> filterMW;              // MW filter

            satYawDataMap satYawData;
            satYawData = satAttitude.getAttitude(gps,satOrbitECI);

            satPCenter.setSatYawData( satYawData );
            windUp.setSatYawData( satYawData );

            gData >> satPCenter             // PCO,PCV of SV
                  >> windUp;                // Wind Up

            gData >> linearAlign            // Phase Code Align
                  >> phaseAlignGPSL1        // Align GPS L1
                  >> phaseAlignGPSL2        // Align GPS L2

                  >> linearPC               // PC
                  >> linearLC               // LC
                  >> prefitForPCE;

            double clock2( Counter::now() );

            satValueMap satClock( basicModel.getSatClockBCE() );


            gData.keepOnlyTypeID( keepTypes );

            currVarSet.clear();

/*
            for(gnssDataMap::iterator gdmIt = gData.begin();
                gdmIt != gData.end();
                ++gdmIt)
            {
                for(sourceDataMap::iterator sdmIt = gdmIt->second.begin();
                    sdmIt != gdmIt->second.end();
                    ++sdmIt)
                {
                    if(sdmIt->first.sourceName != "DJIG") continue;

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

//                    break;
                }

//                break;
            }
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
                    satTypeValueMap stvm( sdmIt->second );

                    Variable var;

                    var = clkStaVar;
                    var.setSource( source );
                    if(currVarSet.find(var) == currVarSet.end())
                        currVarSet.insert( var );

                    var = zwdStaVar;
                    var.setSource( source );
                    if(currVarSet.find(var) == currVarSet.end())
                        currVarSet.insert( var );

                    var = updStaVar;
                    var.setSource( source );
                    wlVarSet.insert( var );

                    for(satTypeValueMap::iterator stvmIt = stvm.begin();
                        stvmIt != stvm.end();
                        ++stvmIt)
                    {
                        SatID sat( stvmIt->first );
                        typeValueMap tvm( stvmIt->second );

                        var = clkSatVar;
                        var.setSatellite( sat );
                        if(currVarSet.find(var) == currVarSet.end())
                            currVarSet.insert( var );

                        var = updSatVar;
                        var.setSatellite( sat );
                        wlVarSet.insert( var );

                        if( tvm.find(TypeID::prefitC12ForPCE) == tvm.end() ||
                            tvm.find(TypeID::prefitL12ForPCE) == tvm.end() )
                            continue;

                        var = ambi12Var;
                        var.setSource( source );
                        var.setSatellite( sat );
                        currVarSet.insert( var );

                        numObs += 2;
                    }
                }
            }

            // set up index
            int id(0);
            for(VariableSet::iterator currIt = currVarSet.begin();
                currIt != currVarSet.end();
                ++currIt)
            {
                Variable& currVar( (Variable&)(*currIt) );

                currVar.setCurrentIndex( id++ );

                VariableSet::iterator prevIt( prevVarSet.find(*currIt) );

                // old variable
                if(prevIt != prevVarSet.end())
                {
                    currVar.setPreviousIndex( prevIt->getCurrentIndex() );

                    TypeID type( currIt->getType() );

                    // ambiguity
                    if(type == TypeID::BLC12)
                    {
                        SourceID source( currIt->getSource() );
                        SatID sat( currIt->getSatellite() );
                        TypeID CSL1( TypeID::CSL1 );

                        // cycle slip
                        if(gData.getValue(gps,source,sat,CSL1) != 0.0)
                        {
                            currVar.setPreviousIndex( -1 );
                        }
                    }
                }
            }


            // state and covar
            int numUnk( currVarSet.size() );

            currState.resize(numUnk);
            for(int i=0; i<numUnk; ++i) currState(i) = 0.0;
            currCovar = Eigen::MatrixXd::Zero(numUnk,numUnk);

            VariableSet tempVarSet( currVarSet );

            for(VariableSet::iterator it1 = currVarSet.begin();
                it1 != currVarSet.end();
                ++it1)
            {
                Variable var1( *it1 );

                int prevIndex1( var1.getPreviousIndex() );
                int currIndex1( var1.getCurrentIndex() );

                TypeID type( var1.getType() );

                // new variable
                if(prevIndex1 == -1)
                {
                    currState(currIndex1) = 0.0;
                    currCovar(currIndex1,currIndex1) = var1.getInitialVariance();

                    if(type == TypeID::dcdtSat)
                    {
                        SatID sat( var1.getSatellite() );

                        if(satClock.find(sat) != satClock.end())
                        {
                            currState(currIndex1) = satClock[sat];
                            currCovar(currIndex1,currIndex1) = 1e+2*1e+2;
                        }
                    }
                }
                // old variable
                else
                {
                    currState(currIndex1) = prevState(prevIndex1);
                    currCovar(currIndex1,currIndex1) = prevCovar(prevIndex1,prevIndex1);
                }

                tempVarSet.erase( *it1 );

                for(VariableSet::iterator it2 = tempVarSet.begin();
                    it2 != tempVarSet.end();
                    ++it2)
                {
                    Variable var2( *it2 );

                    int prevIndex2( var2.getPreviousIndex() );
                    int currIndex2( var2.getCurrentIndex() );

                    // new element
                    if(prevIndex1 == -1 || prevIndex2 == -1)
                    {
                        currCovar(currIndex1,currIndex2) = 0.0;
                        currCovar(currIndex2,currIndex1) = currCovar(currIndex1,currIndex2);
                    }
                    // old element
                    else
                    {
                        currCovar(currIndex1,currIndex2) = prevCovar(prevIndex1,prevIndex2);
                        currCovar(currIndex2,currIndex1) = currCovar(currIndex1,currIndex2);
                    }
                }
            }


            int numWL( wlVarSet.size() );

            wlState.resize(numWL);
            for(int i=0; i<numWL; ++i) wlState(i) = 0.0;
            wlCovar = Eigen::MatrixXd::Zero(numWL,numWL);

            for(VariableSet::iterator it = wlVarSet.begin();
                it != wlVarSet.end();
                ++it)
            {
                Variable var( *it );

                SourceID source( var.getSource() );
                SatID sat( var.getSatellite() );

                int index( var.getCurrentIndex() );
            }

//            for(VariableSet::iterator it = currVarSet.begin();
//                it != currVarSet.end();
//                ++it)
//            {
//                int currIndex(-1), prevIndex(-1);
//                currIndex = it->getCurrentIndex();
//                prevIndex = it->getPreviousIndex();
//
//                cout << setw(15) << it->getType()
//                     << setw(5) << currIndex
//                     << setw(5) << prevIndex
//                     << setw(20) << currState(currIndex)
//                     << setw(20) << currCovar(currIndex,currIndex)
//                     << endl;
//            }
//
//            break;

            for(VariableSet::iterator it = wlVarSet.begin();
                it != wlVarSet.end();
                ++it)
            {
                int index(-1);
                index = it->getCurrentIndex();

                cout << setw(15) << it->getType()
                     << setw(5) << index
                     << setw(20) << currState(index)
                     << setw(20) << currCovar(index,index)
                     << endl;
            }

            break;

            //////// end of filter initialization ////////


            //////// start of measurement update ////////

            // H
            Eigen::MatrixXd H = Eigen::MatrixXd::Zero(numObs+1,numUnk);

            // omc
            Eigen::VectorXd omc(numObs+1);
            for(int i=0; i<numObs+1; ++i) omc(i) = 0.0;

            // R
            Eigen::MatrixXd R = Eigen::MatrixXd::Zero(numObs+1,numObs+1);


            int row(0);

            double obs(0.0), com(0.0);
            for(SatIDSet::iterator it = satSet.begin();
                it != satSet.end();
                ++it)
            {
                SatID sat(*it);

                if(sat.system != SatID::systemGPS) continue;

                Variable var;
                var = clkSatVar;
                var.setSatellite( sat );

                VariableSet::iterator currIt( currVarSet.find(var) );
                if(currIt == currVarSet.end()) continue;

                int index( currIt->getCurrentIndex() );

                if(satClock.find(sat) == satClock.end()) continue;

                double clock( satClock[sat] );

                obs += clock;
                com += currState(index);
                H(row,index) = 1.0;
            }

            omc(row) = obs - com;
            R(row,row) = 1e+1;

            row++;


            double a12( GAMMA_GPS_L1L2/(GAMMA_GPS_L1L2 - 1.0) );
            double b12( -1.0          /(GAMMA_GPS_L1L2 - 1.0) );
            double X( a12*a12 + b12*b12 );

            double cdtr(0.0), zwdr(0.0);
            double cdts(0.0);
            double ambi(0.0);
            double prefitC(0.0), prefitL(0.0);

            for(gnssDataMap::iterator gdmIt = gData.begin();
                gdmIt != gData.end();
                ++gdmIt)
            {
                for(sourceDataMap::iterator sdmIt = gdmIt->second.begin();
                    sdmIt != gdmIt->second.end();
                    ++sdmIt)
                {
                    SourceID source( sdmIt->first );

                    Variable var;

                    int clkStaID(-1);
                    var = clkStaVar;
                    var.setSource( source );
                    if(currVarSet.find(var) == currVarSet.end()) continue;
                    clkStaID = currVarSet.find(var)->getCurrentIndex();
                    cdtr = currState(clkStaID);

                    int zwdStaID(-1);
                    var = zwdStaVar;
                    var.setSource( source );
                    if(currVarSet.find(var) == currVarSet.end()) continue;
                    zwdStaID = currVarSet.find(var)->getCurrentIndex();
                    zwdr = currState(zwdStaID);

                    double clock( sourceClock[source] );

                    for(satTypeValueMap::iterator stvmIt = sdmIt->second.begin();
                        stvmIt != sdmIt->second.end();
                        ++stvmIt)
                    {
                        SatID sat( stvmIt->first );
                        typeValueMap tvm( stvmIt->second );

                        if(sat.system != SatID::systemGPS) continue;

                        if(tvm.find(TypeID::wetMap) == tvm.end()) continue;
                        double wmf( tvm[TypeID::wetMap] );

                        if(tvm.find(TypeID::weight) == tvm.end()) continue;
                        double weightC( tvm[TypeID::weight]*1e+0 );
                        double weightL( tvm[TypeID::weight]*1e+4 );

                        int clkSatID(-1);
                        var = clkSatVar;
                        var.setSatellite( sat );
                        if(currVarSet.find(var) == currVarSet.end()) continue;
                        clkSatID = currVarSet.find(var)->getCurrentIndex();

                        int ambi12ID(-1);
                        var = ambi12Var;
                        var.setSource( source );
                        var.setSatellite( sat );
                        if(currVarSet.find(var) == currVarSet.end()) continue;
                        ambi12ID = currVarSet.find(var)->getCurrentIndex();

                        if( tvm.find(TypeID::prefitC12ForPCE) == tvm.end() ||
                            tvm.find(TypeID::prefitL12ForPCE) == tvm.end() )
                            continue;

                        prefitC = tvm[TypeID::prefitC12ForPCE];
                        prefitL = tvm[TypeID::prefitL12ForPCE];

                        cdts = currState(clkSatID);
                        ambi = currState(ambi12ID);

                        H(row+0,clkStaID) = +1.0;
                        H(row+0,zwdStaID) = +wmf;
                        H(row+0,clkSatID) = -1.0;

                        H(row+1,clkStaID) = +1.0;
                        H(row+1,zwdStaID) = +wmf;
                        H(row+1,clkSatID) = -1.0;
                        H(row+1,ambi12ID) = -1.0;

                        omc(row+0) = prefitC - clock - cdtr - wmf*zwdr + cdts;
                        omc(row+1) = prefitL - clock - cdtr - wmf*zwdr + cdts + ambi;

                        R(row+0,row+0) = X/weightC*9e-2;
                        R(row+1,row+1) = X/weightL*9e-2;

                        row += 2;

                    } // End of for(satTypeValueMap::iterator stvmIt = ...)

                } // End of for(sourceDataMap::iterator sdmIt = ...)

            } // End of for(gnssDataMap::iterator gdmIt = ...)


            Eigen::MatrixXd D( currCovar*H.transpose() );
            Eigen::MatrixXd S( H*D + R );
            Eigen::MatrixXd U( S.llt().matrixU() );
            Eigen::MatrixXd L( U.transpose() );
            Eigen::MatrixXd E( L.triangularView<Eigen::Lower>().solve(D.transpose()).transpose() );
            Eigen::MatrixXd K( U.triangularView<Eigen::Upper>().solve(E.transpose()).transpose() );

            currState = currState + K*omc;
            currCovar = currCovar - E*E.transpose();

            double clock3( Counter::now() );


            outFile << setw(4) << CivilTime(gps).year
                    << setw(3) << CivilTime(gps).month
                    << setw(3) << CivilTime(gps).day
                    << setw(3) << CivilTime(gps).hour
                    << setw(3) << CivilTime(gps).minute
                    << setprecision(2)
                    << setw(6) << CivilTime(gps).second;

            outFile << setw(4) << numSource
                    << setw(4) << numSat
                    << setw(6) << numUnk
                    << setw(6) << numObs;

            outFile << setw(6) << clock2-clock1
                    << setw(6) << clock3-clock2;

            double clk_com(0.0), clk_ref(0.0);

            for(SatIDSet::iterator it = allSatSet.begin();
                it != allSatSet.end();
                ++it)
            {
                SatID sat( *it );

                int num(0);

                if(satSet.find(sat) != satSet.end())
                {
                    try
                    {
                        Variable var;
                        var = clkSatVar;
                        var.setSatellite( sat );

                        int index(currVarSet.find(var)->getCurrentIndex());

                        num = satSourceNumMap[sat];

                        clk_com = currState(index)/C_MPS;
                        clk_ref = sp3Store.getXvt(sat,gps).clkbias;
                    }
                    catch(...)
                    {
                        clk_com = 0.0;
                        clk_ref = 0.0;
                    }
                }
                else
                {
                    clk_com = 0.0;
                    clk_ref = 0.0;
                }

                outFile << setprecision(3)
                        << setw( 3) << sat.id;

                outFile << setw( 3) << num;

                outFile << setw(15) << clk_com*1e+9;
            }

            outFile << endl;

            //////// end of measment update ////////


            //////// start of time update ////////

            for(VariableSet::iterator it = currVarSet.begin();
                it != currVarSet.end();
                ++it)
            {
                Variable var(*it);

                TypeID type( var.getType() );
                int index( var.getCurrentIndex() );

                if(type == TypeID::dcdtSta)
                {
//                    currState(index) = 0.0;
                    currCovar(index,index) += 1e+2*1e+2;
                }
                else if(type == TypeID::wetMap)
                {
                    currCovar(index,index) += 1e-8*dt;
                }
                else if(type == TypeID::dcdtSat)
                {
//                    currState(index) = 0.0;
                    currCovar(index,index) += 1e+2*1e+2;
                }
            }

            //////// end of time update ////////

            prevVarSet = currVarSet;
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

//            if((gps-gps1) > 120*30.0) break;

            if(gps >= gps2) break;

        }
        catch(...)
        {
            cerr << "processing error." << endl;
            break;
        }

    } // End of 'while( true )'

    outFile.close();

    double clock_end( Counter::now() );

    cerr << fixed;
    cerr << setw(20) << "time elapsed: "
         << setw(10) << setprecision(3)
         << clock_end - clock_start << endl;

    return 0;

} // main()
