#pragma ident "$Id$"

#include <algorithm>

#include "ConfDataReader.hpp"

#include "EOPDataStore2.hpp"
#include "LeapSecStore.hpp"

#include "ReferenceSystem.hpp"

#include "SolarSystem.hpp"

#include "GNSSOrbit.hpp"

#include "RKF78Integrator.hpp"

#include "MSCStore.hpp"

#include "Rinex3NavStream.hpp"
#include "Rinex3NavHeader.hpp"
#include "Rinex3NavData.hpp"

#include "Rinex3EphemerisStore2.hpp"

#include "Rinex3ObsStream.hpp"
#include "Rinex3ObsHeader.hpp"
#include "Rinex3ObsData.hpp"

#include "SP3EphemerisStore.hpp"

#include "NetworkObsStreams.hpp"

#include "DataStructures.hpp"

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


int main(int argc, char* argv[])
{

    double dt(30.0);

    // conf file
    string confFileName("orbit_clock.conf");

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

    confReader.setFallback2Default(true);


    // initial time
    string time;
    try
    {
        time = confReader.getValue("initialTime", "DEFAULT");
    }
    catch(...)
    {
        cerr << "initial time get error." << endl;
        exit(-1);
    }

    int year = asInt( time.substr( 0,4) );
    int mon  = asInt( time.substr( 5,2) );
    int day  = asInt( time.substr( 8,2) );
    int hour = asInt( time.substr(11,2) );
    int min  = asInt( time.substr(14,2) );
    int sec  = asDouble( time.substr(17,2) );

    CivilTime civilTime(year,mon,day,hour,min,sec, TimeSystem::GPS);
    CommonTime gps0( civilTime.convertToCommonTime() );

//    cout << "initial time: " << civilTime << endl;


    // eop file
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

    EOPDataStore2 eopStore;

    try
    {
        eopStore.loadIERSFile(eopFileName);
    }
    catch(...)
    {
        cerr << "eop file '" << eopFileName
             << "' load error." << endl;
        exit(-1);
    }


    // leap second file
    string lsFileName;

    try
    {
        lsFileName = confReader.getValue("lsFileName", "DEFAULT");
    }
    catch(...)
    {
        cerr << "leap second file name get error." << endl;
        exit(-1);
    }

    LeapSecStore lsStore;

    try
    {
        lsStore.loadFile(lsFileName);
    }
    catch(...)
    {
        cerr << "leap second file '" << lsFileName
             << "' load error." << endl;
        exit(-1);
    }


    // reference system
    ReferenceSystem refSys;
    refSys.setEOPDataStore(eopStore);
    refSys.setLeapSecStore(lsStore);


    // sp3 file
    string sp3FileListName;

    try
    {
        sp3FileListName = confReader.getValue("sp3FileListName", "DEFAULT");
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


    // clk file
    string clkFileListName;

    try
    {
        clkFileListName = confReader.getValue("clkFileListName", "DEFAULT");
    }

    catch(...)
    {
        cerr << "clk file list name get error." << endl;
        exit(-1);
    }

    vector<string> clkFileListVec;

    ifstream clkFileListStream;

    clkFileListStream.open(clkFileListName.c_str());

    if(!clkFileListStream)
    {
        cerr << "clk file list '" << clkFileListName
             << "' open error." << endl;
        exit(-1);
    }

    string clkFile;
    while(clkFileListStream >> clkFile)
    {
        clkFileListVec.push_back(clkFile);
    }

    clkFileListStream.close();

    if(clkFileListVec.size() == 0)
    {
        cerr << "clk file list is empty." << endl;
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

    vector<string>::const_iterator clk_it = clkFileListVec.begin();
    while(clk_it != clkFileListVec.end())
    {
        clkFile = (*clk_it);

        try
        {
            sp3Store.loadRinexClockFile(clkFile);
        }
        catch(...)
        {
            cerr << "clk file '" << clkFile
                 << "' load error." << endl;

            ++clk_it;
            continue;
        }

        ++clk_it;
    }


    // nav file
    string navFileListName;

    try
    {
        navFileListName = confReader.getValue("navFileListName", "DEFAULT");
    }
    catch(...)
    {
        cerr << "nav file list name get error." << endl;
        exit(-1);
    }

    vector<string> navFileListVec;

    ifstream navFileListStream;

    navFileListStream.open(navFileListName.c_str());

    if(!navFileListStream)
    {
        cerr << "nav file list '" << navFileListName
             << "' open error." << endl;
        exit(-1);
    }

    string navFile;
    while(navFileListStream >> navFile)
    {
        navFileListVec.push_back(navFile);
    }

    navFileListStream.close();

    if(navFileListVec.size() == 0)
    {
        cerr << "nav file list is empty." << endl;
        exit(-1);
    }


    Rinex3EphemerisStore2 bceStore;

    vector<string>::const_iterator nav_it = navFileListVec.begin();
    while(nav_it != navFileListVec.end())
    {
        string navFile = (*nav_it);

        try
        {
            bceStore.loadFile(navFile);
        }
        catch(...)
        {
            cerr << "nav file '" << navFile
                 << "' load error." << endl;

            ++nav_it;

            continue;
        }

        ++nav_it;
    }


    Matrix<double> t2cRaw(3,3,0.0), t2cDot(3,3,0.0);


    // initial conditions
    CommonTime utc0( refSys.GPS2UTC(gps0) );
    CommonTime tt0( refSys.GPS2TT(gps0) );

    // orbit
    satVectorMap satOrbit0;

    // srpc
    satVectorMap satSRPC0;

    // clock
    satValueMap satClock0;

    t2cRaw = refSys.T2CMatrix(utc0);
    t2cDot = refSys.dT2CMatrix(utc0);

    for(int i=1; i<=MAX_PRN_GPS; ++i)
    {
        SatID sat(i,SatID::systemGPS);

        Vector<double> rsat_t(3,0.0), vsat_t(3,0.0);
        Vector<double> rsat_c(3,0.0), vsat_c(3,0.0);

        Vector<double> srpc0(5,0.0);
        srpc0(0) = -1.0;

        double clock0(0.0);

        try
        {
            rsat_t = sp3Store.getXvt(sat,gps0).x.toVector();
            vsat_t = sp3Store.getXvt(sat,gps0).v.toVector();
            clock0 = sp3Store.getXvt(sat,gps0).clkbias * C_MPS;
        }
        catch(...)
        {
            cerr << "initial conditions of " << sat << " get error." << endl;
            continue;
        }

        rsat_c = t2cRaw * rsat_t;
        vsat_c = t2cRaw * vsat_t + t2cDot * rsat_t;

        // (r, v, dr/dr0, dr/dv0, dv/dr0, dv/dv0, dr/dp0, dv/dp0)
        Vector<double> orbit0(72,0.0);
        orbit0( 6) = 1.0; orbit0(10) = 1.0; orbit0(14) = 1.0; // dr/dr0
        orbit0(33) = 1.0; orbit0(37) = 1.0; orbit0(41) = 1.0; // dv/dv0

        orbit0(0) = rsat_c[0]; orbit0(1) = rsat_c[1]; orbit0(2) = rsat_c[2];
        orbit0(3) = vsat_c[0]; orbit0(4) = vsat_c[1]; orbit0(5) = vsat_c[2];

        satOrbit0[sat] = orbit0;
        satSRPC0[sat] = srpc0;
        satClock0[sat] = clock0;
    }


    // jpl ephemeris file
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


    // solar system
    SolarSystem solSys;

    try
    {
        solSys.initializeWithBinaryFile(jplFileName);
    }
    catch(...)
    {
        cerr << "solar system initialize error." << endl;
        exit(-1);
    }


    // egm file
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


    // earth gravitation
    EGM08Model egm(12,12);

    try
    {
        egm.loadFile(egmFileName);
    }
    catch(...)
    {
        cerr << "egm file '" << egmFileName << "' load error." << endl;
        exit(-1);
    }

    egm.setReferenceSystem(refSys);


    // solid tide
    EarthSolidTide solidTide;
    solidTide.setReferenceSystem(refSys);
    solidTide.setSolarSystem(solSys);

    egm.setEarthSolidTide(solidTide);


    // ocean tide file
    string otFileName;

    try
    {
        otFileName = confReader.getValue("otFileName", "DEFAULT");
    }
    catch(...)
    {
        cerr << "ocean tide file name get error." << endl;
        exit(-1);
    }


    // ocean tide
    EarthOceanTide oceanTide(8,8);
    oceanTide.setReferenceSystem(refSys);

    try
    {
        oceanTide.loadFile(otFileName);
    }
    catch(...)
    {
        cerr << "ocean tide file '" << otFileName << "' load error." << endl;
        exit(-1);
    }

    egm.setEarthOceanTide(oceanTide);


    // pole tide
    EarthPoleTide poleTide;
    poleTide.setReferenceSystem(refSys);

    egm.setEarthPoleTide(poleTide);


    // third body gravitation
    ThirdBody thd;
    thd.setReferenceSystem(refSys);
    thd.setSolarSystem(solSys);
    thd.enableAllPlanets();


    // solar radiation pressure
    ECOM1Model srp;
    srp.setReferenceSystem(refSys);
    srp.setSolarSystem(solSys);


    // relativity
    Relativity rel;


    // gnss orbit
    GNSSOrbit gnss;
    gnss.setEGMModel( egm );
    gnss.setThirdBody( thd );
    gnss.setSRPModel( srp );
    gnss.setRelativity( rel );


    // rkf78 integrator
    RKF78Integrator rkf78;
    rkf78.setStepSize( dt );


    // msc file
    string mscFileName;
    try
    {
        mscFileName = confReader.getValue("mscFileName", "DEFAULT");
    }
    catch(...)
    {
        cerr << "get msc file name error." << endl;
        exit(-1);
    }

    MSCStore mscStore;

    try
    {
        mscStore.loadFile( mscFileName.c_str() );
    }
    catch(...)
    {
        cerr << "msc file '" << mscFileName << "' open error." << endl;
        exit(-1);
    }


    // dcb file
    string dcbFileName;

    try
    {
        dcbFileName = confReader.getValue("dcbFileName", "DEFAULT");
    }
    catch(...)
    {
        cerr << "dcb file name get error." << endl;
        exit(-1);
    }

    DCBDataReader dcbReader;

    try
    {
        dcbReader.open( dcbFileName );
    }
    catch(...)
    {
        cerr << "dcb file '" << dcbFileName << "' open error." << endl;
        exit(-1);
    }


    // atx file
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

    AntexReader antexReader;

    try
    {
        antexReader.open( atxFileName );
    }
    catch(...)
    {
        cerr << "atx file '" << atxFileName << "' open error." << endl;
        exit(-1);
    }


    // blq file
    BLQDataReader blqStore;

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
        blqStore.open( blqFileName );
    }
    catch(...)
    {
        cerr << "blq file '" << blqFileName << "' open error." << endl;
        exit(-1);
    }


    // gpt2 file
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


    // obs file
    NetworkObsStreams obsStreams;

    string obsFileListName;
    try
    {
        obsFileListName = confReader.getValue("obsFileListName", "DEFAULT");
    }
    catch(...)
    {
        cerr << "obs file list name get error." << endl;
        exit(-1);
    }

    ifstream obsFileListStream( obsFileListName.c_str() );

    if( !obsFileListStream.is_open() )
    {
        cerr << "obs file list '" << obsFileListName << "' open error." << endl;
        exit(-1);
    }


    double clock_start( Counter::now() );


    SourceIDSet allSourceSet;

    map<SourceID, Triple> sourceMonumentMap;
    map<SourceID, Antenna> sourceAntennaMap;


    // ObsID on GPS L1 Frequency
    RinexObsID GC1C("GC1C"), GC1W("GC1W"), GL1C("GL1C");
    // ObsID on GPS L2 Frequency
    RinexObsID GC2W("GC2W"), GL2W("GL2W");

    string obsFile;
    while(obsFileListStream >> obsFile)
    {
        Rinex3ObsStream ros;
        ros.exceptions(ios::failbit);

        try
        {
            ros.open(obsFile.c_str(), ios::in);
        }
        catch(...)
        {
            cerr << "obs file '" << obsFile << "' open error." << endl;
            ros.close(); continue;
        }

        Rinex3ObsHeader roh;

        try
        {
            ros >> roh;
        }
        catch(...)
        {
            cerr << "obs header of '" << obsFile << "' read error." << endl;
            ros.close(); continue;
        }


        SourceID source(SourceID::Mixed,roh.markerName,roh.markerNumber);

        string station( source.sourceName.substr(0,4) );
        transform(station.begin(),station.end(),station.begin(),::toupper);


        // check if station is included in MSC file
        gps0.setTimeSystem(TimeSystem::Unknown);

        MSCData mscData;

        try
        {
            mscData = mscStore.findMSC(station, gps0);
        }
        catch(...)
        {
            cerr << "station '" << station << "' not exist in MSC file." << endl;
            ros.close(); continue;
        }

        gps0.setTimeSystem(TimeSystem::GPS);


        // check if current station exist in BLQ file
        if( !blqStore.isValid(station) )
        {
            cerr << "station '" << station << "' not exist in BLQ file." << endl;
            ros.close(); continue;
        }


        // check if required signals exist in OBS file
        map< string,vector<RinexObsID> > mapObsTypes;
        mapObsTypes = roh.mapObsTypes;

        vector<RinexObsID> roi;

        if(mapObsTypes.find("G") != mapObsTypes.end())
        {
            roi = mapObsTypes["G"];

            if( find(roi.begin(),roi.end(),GC1C) == roi.end() &&
                find(roi.begin(),roi.end(),GC1W) == roi.end() )
            {
                ros.close(); continue;
            }

            if( find(roi.begin(),roi.end(),GC2W) == roi.end() ||
                find(roi.begin(),roi.end(),GL1C) == roi.end() ||
                find(roi.begin(),roi.end(),GL2W) == roi.end() )
            {
                ros.close(); continue;
            }
        }

        ros.close();

        // now, we can add OBS file to OBS streams
        if( !obsStreams.addRinexObsFile( obsFile ) )
        {
            cerr << "obs file '" << obsFile << "' add error." << endl;
            continue;
        }

        string antennaModel( roh.antType );
        Antenna antenna;

        try
        {
            antenna = antexReader.getAntenna( antennaModel );
        }
        catch(ObjectNotFound& notFound)
        {
            antennaModel.replace(16,4,"NONE");
            antenna = antexReader.getAntenna( antennaModel );
        }


        allSourceSet.insert(source);
        sourceMonumentMap[source] = roh.antennaDeltaHEN;
        sourceAntennaMap[source] = antenna;
    }

    obsFileListStream.close();


    // CC2NONCC
    CC2NONCC cc2noncc;
    cc2noncc.setDCBDataReader(dcbReader);


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
    ComputeLinear prefitPC;
    prefitPC.addLinear(SatID::systemGPS, linearComb.pc12PrefitOfGPS);

    // ComputeLinear, prefit PC for POD
    ComputeLinear prefitPCForPOD;
    prefitPCForPOD.addLinear(SatID::systemGPS, linearComb.pc12PrefitOfGPSForPOD);

    // ComputeLinear, prefit LC for POD
    ComputeLinear prefitLCForPOD;
    prefitLCForPOD.addLinear(SatID::systemGPS, linearComb.lc12PrefitOfGPSForPOD);


    // MWCSDetector
    MWCSDetector markCSMW;
    markCSMW.setMaxNumLambdas( 3.0 );


    // LICSDetector
    LICSDetector markCSLI;
    markCSLI.setMinThreshold( 0.240 );
    markCSLI.setLIDrift( 0.002 );


    // SatArcMarker
    SatArcMarker markArc;
    markArc.setDeleteUnstableSats(true);
    markArc.setUnstablePeriod(30.0*5+1.0);



    // BasicModel
    BasicModel basicModel;
    basicModel.setEphStore(sp3Store);
    basicModel.setMSCStore(mscStore);
    basicModel.setMinElev(10.0);
    basicModel.setDefaultObs(SatID::systemGPS, TypeID::PC);


    // ComputeElevWeights
    ComputeElevWeights elevWeights;


    // ComputeTropModel
    ViennaTropModel viennaTM;
    viennaTM.loadFile(gpt2FileName);
    ComputeTropModel computeTM;
    computeTM.setTropModel(viennaTM);
    computeTM.setMSCStore(mscStore);


    // ComputeStaTides
    ComputeStaTides staTides;
    staTides.setBLQDataReader(blqStore);
    staTides.setReferenceSystem(refSys);
    staTides.setSolarSystem(solSys);


    // CorrectObservables
    CorrectObservables correctObs;
    correctObs.setMSCStore(mscStore);
    correctObs.setTideCorr(staTides);
    correctObs.setSourceMonument(sourceMonumentMap);
    correctObs.setSourceAntenna(sourceAntennaMap);


    // GravitationalDelay
    GravitationalDelay gravDelay;
    gravDelay.setMSCStore(mscStore);


    // ComputeStaClock
    ComputeStaClock computeStaClock;
    computeStaClock.setMinSatOfGPS(4);
    computeStaClock.setMaxRMSOfGPS(3.0);


    // ComputeSatPCenter
    ComputeSatPCenter satPCenter;
    satPCenter.setReferenceSystem(refSys);
    satPCenter.setSolarSystem(solSys);
    satPCenter.setEphStore(sp3Store);
    satPCenter.setMSCStore(mscStore);
    satPCenter.setAntexReader(antexReader);


    // ComputeWindUp
    ComputeWindUp windUp;
    windUp.setReferenceSystem(refSys);
    windUp.setSolarSystem(solSys);
    windUp.setEphStore(sp3Store);
    windUp.setMSCStore(mscStore);


    // PhaseCodeAlignment
    PhaseCodeAlignment phaseAlignL1;
    phaseAlignL1.setSatSystem(SatID::systemGPS);
    phaseAlignL1.setCodeType(TypeID::Q1);
    phaseAlignL1.setPhaseType(TypeID::L1);
    phaseAlignL1.setPhaseWavelength(L1_WAVELENGTH_GPS);
    phaseAlignL1.setUseSatArc(true);

    PhaseCodeAlignment phaseAlignL2;
    phaseAlignL2.setSatSystem(SatID::systemGPS);
    phaseAlignL2.setCodeType(TypeID::Q2);
    phaseAlignL2.setPhaseType(TypeID::L2);
    phaseAlignL2.setPhaseWavelength(L2_WAVELENGTH_GPS);
    phaseAlignL2.setUseSatArc(true);


    // ComputePartials
    ComputePartials computePartials;
    computePartials.setReferenceSystem(refSys);
    computePartials.setMSCStore(mscStore);


    gnssDataMap gData;

    SourceIDSet sourceSet;
    SatIDSet satSet;

    int numSource(0);
    int numSat(0);

    sourceValueMap sourceClock;
    satValueMap satClock;

    sourceVectorMap sourcePosECI;

    VariableVector prevVarVec, currVarVec;

    Vector<double> prevState, currState;
    Matrix<double> prevCovar, currCovar;


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
    Variable clkSat(TypeID::dcdtSat);
    clkSat.setSourceIndexed( false );
    clkSat.setSatIndexed( true );
    clkSat.setInitialVariance( 1e+2*1e+2 );

    VariableVector svVarVec;
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

    keepTypes.insert(TypeID::prefitCForPOD);
    keepTypes.insert(TypeID::prefitLForPOD);
    keepTypes.insert(TypeID::weight);
    keepTypes.insert(TypeID::wetMap);
    keepTypes.insert(TypeID::CSL1);
    keepTypes.insert(TypeID::rho);


    SatIDSet allSatSet;

    for(int i=1; i<=MAX_PRN_GPS; ++i)
    {
        if(i == 4) continue;

        SatID sat(i,SatID::systemGPS);
        allSatSet.insert(sat);
    }


    CommonTime gps, utc, tt;

//    ofstream fres("residual.txt");
//    fres << fixed << setprecision(3);

    cout << fixed;

    // process epoch by epoch
    while( true )
    {

        obsStreams.readEpochData(gps, gData);

        gData.keepOnlySatID(allSatSet);

        gps = gData.begin()->first;

//        if((gps-gps0) >= 600 * 30.0) break;

        try
        {
            double clock1( Counter::now() );


            utc = refSys.GPS2UTC( gps );
            tt = refSys.GPS2TT( gps );

            t2cRaw = refSys.T2CMatrix(utc);
            t2cDot = refSys.dT2CMatrix(utc);


            SatIDSet satRejectedSet;

            satVectorMap satPosECI, satVelECI;

            for(SatIDSet::iterator it = allSatSet.begin();
                it != allSatSet.end();
                ++it)
            {
                SatID sat( *it );

                Vector<double> rsat_t(3,0.0), vsat_t(3,0.0);
                Vector<double> rsat_c(3,0.0), vsat_c(3,0.0);

                try
                {
                    rsat_t = sp3Store.getXvt(sat,gps).x.toVector();
                    vsat_t = sp3Store.getXvt(sat,gps).v.toVector();

                    rsat_c = t2cRaw * rsat_t;
                    vsat_c = t2cRaw * vsat_t + t2cDot * rsat_t;

                    satPosECI[sat] = rsat_c;
                    satVelECI[sat] = vsat_c;
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
                  >> basicModel;            // r,v,cdt of SV
                                            // rho,elev,azim of RX-SV

            satClock = basicModel.getSatClockMap();

            gData >> elevWeights            // weight
                  >> computeTM              // troposphere
                  >> correctObs             // disp,ARP,PCO,PCV of RX
                  >> satPCenter             // PCO,PCV of SV
                  >> gravDelay              // gravitational delay
                  >> linearPC
                  >> prefitPC
                  >> computeStaClock;

            sourceClock = computeStaClock.getSourceClockMap();

            computePartials.setSourceClockMap( sourceClock );

            gData >> linearCS               // Cycle Slip
                  >> markCSMW               // MW detector
                  >> markCSLI               // LI detector
                  >> markArc                // obs arc

                  >> computePartials        // partials
                  >> windUp                 // Wind Up

                  >> linearAlign            // Phase Code Align
                  >> phaseAlignL1           // Align L1
                  >> phaseAlignL2           // Align L2

                  >> linearPC               // PC
                  >> linearLC               // LC
                  >> prefitPCForPOD
                  >> prefitLCForPOD;

            sourcePosECI = computePartials.getSourcePosECIMap();

            double clock2( Counter::now() );


            gData.keepOnlyTypeID( keepTypes );

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

                        if(tvm.find(TypeID::prefitLForPOD) == tvm.end()) continue;

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
            int numUnknown( currVarVec.size() );
            currState.resize(numUnknown,0.0);
            currCovar.resize(numUnknown,numUnknown,0.0);

            for(int i=0; i<numUnknown; ++i)
            {
                Variable var1( currVarVec[i] );

                int prevIndex1( var1.getPreviousIndex() );
                int currIndex1( var1.getCurrentIndex() );

                TypeID type( var1.getType() );

                // new variable
                if(prevIndex1 == -1)
                {
                    if(type == TypeID::dcdtSta)
                    {
                        currState(currIndex1) = 0.0;
                    }
                    else if(type == TypeID::wetMap)
                    {
                        currState(currIndex1) = 0.0;
                    }
                    else if(type == TypeID::dcdtSat)
                    {
                        SatID sat( var1.getSatellite() );
                        currState(currIndex1) = satClock[sat];
                    }
                    else if(type == TypeID::BLC)
                    {
                        currState(currIndex1) = 0.0;
                    }

                    currCovar(currIndex1,currIndex1) = var1.getInitialVariance();
                }
                // old variable
                else
                {
                    currState(currIndex1) = prevState(prevIndex1);
                    currCovar(currIndex1,currIndex1) = prevCovar(prevIndex1,prevIndex1);
                }

                for(int j=i+1; j<numUnknown; ++j)
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

            //////// end of filter initialization ////////


            //////// start of clock constraint ////////

            double obs(0.0), com(0.0);

            // omc, weight
            double omc(0.0), weight(0.0);

            // h
            Vector<double> h(numUnknown,0.0);

            for(SatIDSet::iterator it = satSet.begin();
                it != satSet.end();
                ++it)
            {
                SatID sat(*it);

                if(satIndexMap.find(sat) == satIndexMap.end()) continue;

                int index( satIndexMap[sat] );

                if(satClock.find(sat) == satClock.end()) continue;

                obs += satClock[sat];
                com += currState(index+0);
                h(index+0) = 1.0;
            }

            omc = obs - com;

            // p * h'
            Vector<double> pht(numUnknown,0.0);
            for(int i=0; i<numUnknown; ++i)
            {
                for(int j=0; j<numUnknown; ++j)
                {
                    if(h(j) != 0.0) pht(i) += currCovar(i,j) * h(j);
                }
            }

            // h * p * h'
            double hpht(0.0);
            for(int i=0; i<numUnknown; ++i)
            {
                if(h(i) != 0.0) hpht += h(i) * pht(i);
            }

            weight = 2.0;

            // beta
            double beta(0.0);
            beta = 1.0/weight + hpht;

            // kalman gain
            Vector<double> gamma(numUnknown,0.0);
            gamma = pht/beta;

            // state update
            currState = currState + gamma*omc;

            // covariance update
#pragma omp parallel for
            for(int i=0; i<numUnknown; ++i)
            {
                currCovar(i,i) = currCovar(i,i) - gamma(i)*pht(i);

                for(int j=i+1; j<numUnknown; ++j)
                {
                    currCovar(i,j) = currCovar(i,j) - gamma(i)*pht(j);
                    currCovar(j,i) = currCovar(i,j);
                }
            }

            //////// end of clock constraint ////////


            //////// start of measurement update ////////

            double cdtSta(0.0), zwdSta(0.0);
            Vector<double> rSat(3,0.0), vSat(3,0.0);
            double cdtSat(0.0);
            double ambi(0.0);
            double prefit(0.0), postfit(0.0);

            double r(0.0), r2(0.0), r3(0.0), k(0.0);
            double rho(0.0);

            Vector<double> posSatECI(3,0.0);

//            fres << CivilTime(gps) << endl;

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

                        int svIndex( satIndexMap[sat] );
                        int rxsvIndex( sourceSatIndexMap[source][sat] );

                        if(tvm.find(TypeID::wetMap) == tvm.end()) continue;
                        double wmf( tvm[TypeID::wetMap] );

//                        double dtSta( clock/C_MPS );
//
//                        if(tvm.find(TypeID::rho) == tvm.end()) continue;
//                        double dtRho( tvm[TypeID::rho]/C_MPS );
//
//                        double dT( dtSta + dtRho );
//                        double dT2( dT * dT );
//
//                        rSat = satPosECI[sat];
//                        vSat = satVelECI[sat];
//
//                        r = norm(rSat);
//                        r2 = r * r;
//                        r3 = r2 * r;
//
//                        k = GM_EARTH/r3;
//
//                        posSatECI = rSat - vSat*dT - 0.5*k*rSat*dT2;
//
//                        rho = norm(posSatECI-posSourceECI);

                        if(tvm.find(TypeID::rho) == tvm.end()) continue;
                        double rho( tvm[TypeID::rho] );

//                        fres << setw(5) << source.sourceName.substr(0,4)
//                             << setw(5) << sat
//                             << setw(10) << dT
//                             << setw(15) << tvm[TypeID::rho]
//                             << setw(15) << rho
//                             << endl;

                        //////// code ////////

                        cdtSta = currState(rxIndex+0);
                        zwdSta = currState(rxIndex+1);
                        cdtSat = currState(svIndex+0);

                        // omc, weight, wmf
                        if(tvm.find(TypeID::prefitCForPOD) == tvm.end()) continue;
                        prefit = tvm[TypeID::prefitCForPOD];

                        if(tvm.find(TypeID::weight) == tvm.end()) continue;
                        weight = tvm[TypeID::weight]*1e+0;

                        omc = prefit - rho - cdtSta - wmf*zwdSta + cdtSat;

                        // h
                        h.resize(numUnknown,0.0);
                        h(rxIndex+0) = +1.0;
                        h(rxIndex+1) =  wmf;
                        h(svIndex+0) = -1.0;

                        // p * h'
                        pht.resize(numUnknown,0.0);
                        for(int i=0; i<numUnknown; ++i)
                        {
                            for(int j=0; j<numUnknown; ++j)
                            {
                                if(h(j) != 0.0) pht(i) += currCovar(i,j) * h(j);
                            }
                        }

                        // h * p * h'
                        hpht = 0.0;
                        for(int i=0; i<numUnknown; ++i)
                        {
                            if(h(i) != 0.0) hpht += h(i) * pht(i);
                        }

                        beta = 1.0/weight + hpht;

                        // kalman gain
                        gamma = pht/beta;

                        // state update
                        currState = currState + gamma*omc;

                        // covariance update
#pragma omp parallel for
                        for(int i=0; i<numUnknown; ++i)
                        {
                            currCovar(i,i) = currCovar(i,i) - gamma(i)*pht(i);

                            for(int j=i+1; j<numUnknown; ++j)
                            {
                                currCovar(i,j) = currCovar(i,j) - gamma(i)*pht(j);
                                currCovar(j,i) = currCovar(i,j);
                            }
                        }


                        //////// phase ////////

                        cdtSta = currState(rxIndex+0);
                        zwdSta = currState(rxIndex+1);
                        cdtSat = currState(svIndex+0);
                        ambi = currState(rxsvIndex+0);

                        // omc, weight, wmf
                        if(tvm.find(TypeID::prefitLForPOD) == tvm.end()) continue;
                        prefit = tvm[TypeID::prefitLForPOD];

                        if(tvm.find(TypeID::weight) == tvm.end()) continue;
                        weight = tvm[TypeID::weight]*1e+4;

                        omc = prefit - rho - cdtSta - wmf*zwdSta + cdtSat + ambi;

                        // h
                        h.resize(numUnknown,0.0);
                        h(rxIndex+0) = +1.0;
                        h(rxIndex+1) =  wmf;
                        h(svIndex+0) = -1.0;
                        h(rxsvIndex+0) = -1.0;

                        // p * h'
                        pht.resize(numUnknown,0.0);
                        for(int i=0; i<numUnknown; ++i)
                        {
                            for(int j=0; j<numUnknown; ++j)
                            {
                                if(h(j) != 0.0) pht(i) += currCovar(i,j) * h(j);
                            }
                        }

                        // h * p * h'
                        hpht = 0.0;
                        for(int i=0; i<numUnknown; ++i)
                        {
                            if(h(i) != 0.0) hpht += h(i) * pht(i);
                        }

                        beta = 1.0/weight + hpht;

                        // kalman gain
                        gamma = pht/beta;

                        // state update
                        currState = currState + gamma*omc;

                        // covariance update
#pragma omp parallel for
                        for(int i=0; i<numUnknown; ++i)
                        {
                            currCovar(i,i) = currCovar(i,i) - gamma(i)*pht(i);

                            for(int j=i+1; j<numUnknown; ++j)
                            {
                                currCovar(i,j) = currCovar(i,j) - gamma(i)*pht(j);
                                currCovar(j,i) = currCovar(i,j);
                            }
                        }

                    } // End of for(satTypeValueMap::iterator stvmIt = ...)

                } // End of for(sourceDataMap::iterator sdmIt = ...)

            } // End of for(gnssDataMap::iterator gdmIt = ...)


            // postfit
//            for(gnssDataMap::iterator gdmIt = gData.begin();
//                gdmIt != gData.end();
//                ++gdmIt)
//            {
//                for(sourceDataMap::iterator sdmIt = gdmIt->second.begin();
//                    sdmIt != gdmIt->second.end();
//                    ++sdmIt)
//                {
//                    SourceID source( sdmIt->first );
//
//                    int rxIndex( sourceIndexMap[source] );
//
//                    if(sourcePosECI.find(source) == sourcePosECI.end()) continue;
//                    Vector<double> posSourceECI( sourcePosECI[source] );
//
//                    for(satTypeValueMap::iterator stvmIt = sdmIt->second.begin();
//                        stvmIt != sdmIt->second.end();
//                        ++stvmIt)
//                    {
//                        SatID sat( stvmIt->first );
//
//                        typeValueMap tvm( stvmIt->second );
//
//                        int svIndex( satIndexMap[sat] );
//
//                        cdtSta = currState(rxIndex+0);
//                        zwdSta = currState(rxIndex+1);
//
//                        // sat pos at receive time
//                        rSat(0) = currState(svIndex+0);
//                        rSat(1) = currState(svIndex+1);
//                        rSat(2) = currState(svIndex+2);
//
//                        // sat vel at receive time
//                        vSat(0) = currState(svIndex+3);
//                        vSat(1) = currState(svIndex+4);
//                        vSat(2) = currState(svIndex+5);
//
//                        double dtSta( cdtSta/C_MPS );
//
//                        if(tvm.find(TypeID::rho) == tvm.end()) continue;
//                        double dtRho( tvm[TypeID::rho]/C_MPS );
//
//                        double dT( dtSta + dtRho );
//                        double dT2( dT * dT );
//
//                        r = norm(rSat);
//                        r2 = r * r;
//                        r3 = r2 * r;
//
//                        k = GM_EARTH/r3;
//
//                        // sat pos at transmit time
//                        posSatECI = rSat - vSat*dT - 0.5*k*rSat*dT2;
//
//                        rho = norm(posSatECI-posSourceECI);
//
//                        cdtSat = currState(svIndex+11);
//
//                        double wmf = stvmIt->second[TypeID::wetMap];
//
//                        prefit = stvmIt->second[TypeID::prefitCForPOD];
//
//                        postfit = prefit - rho - cdtSta + cdtSat - wmf*zwdSta;
//                    }
//                }
//            }


            cout << setw(4) << CivilTime(gps).year
                 << setw(4) << CivilTime(gps).month
                 << setw(4) << CivilTime(gps).day
                 << setw(4) << CivilTime(gps).hour
                 << setw(4) << CivilTime(gps).minute
                 << setprecision(3)
                 << setw(10) << CivilTime(gps).second;

            cout << setw(5) << numSource
                 << setw(5) << numSat
                 << setw(5) << numUnknown;

            Vector<double> rcom_c(3,0.0), rcom_t(3,0.0);
            Vector<double> rref_t(3,0.0), rref_c(3,0.0);
            Vector<double> vref_t(3,0.0), vref_c(3,0.0);

            for(SatIDSet::iterator it = allSatSet.begin();
                it != allSatSet.end();
                ++it)
            {
                SatID sat( *it );

                int num(0);
                double com(0.0), ref(0.0);

                if(satSet.find(sat) != satSet.end())
                {
                    try
                    {
                        int index(satIndexMap[sat]);

                        num = satSourceNumMap[sat];

                        com = currState(index+0)/C_MPS;
                        ref = sp3Store.getXvt(sat,gps).clkbias;
                    }
                    catch(...)
                    {
                        num = 0;
                        com = 0.0;
                        ref = 0.0;
                    }
                }

                cout << setprecision(3)
                     << setw(5) << sat.id
                     << setw(5) << num
                     << setw(15) << com * 1e9
                     << setw(15) << ref * 1e9;
            }
            cout << endl;

            double clock4( Counter::now() );

            //////// end of measment update ////////


            //////// start of time update ////////

            double clock5( Counter::now() );


            // source-related variable
            for(SourceIDSet::iterator it = sourceSet.begin();
                it != sourceSet.end();
                ++it)
            {
                SourceID source( *it );

                int index(sourceIndexMap[source]);

//                currState(index+0) = 0.0;

                currCovar(index+0,index+0) += 1e+2*dt;
                currCovar(index+1,index+1) += 3e-9*dt;
            }


            // sat-related variable
            for(SatIDSet::iterator it = satSet.begin();
                it != satSet.end();
                ++it)
            {
                SatID sat( *it );

                int index(satIndexMap[sat]);

                currCovar(index+0,index+0) += 1e-1*dt;
            }

            for(int i=0; i<numUnknown; ++i)
            {
                for(int j=i+1; j<numUnknown; ++j)
                {
                    currCovar(j,i) = currCovar(i,j);
                }
            }

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
        }
        catch(...)
        {
            cerr << "processing error." << endl;
            break;
        }

    } // End of 'while( obsStreams.readEpochData(gData) )'

    double clock_end( Counter::now() );

    cerr << setw(20) << "time elapsed: "
         << setw(10) << setprecision(3)
         << clock_end - clock_start << endl;

//    fres.close();

    return 0;

}
