#pragma ident "$Id$"

#include <algorithm>

#include "ConfDataReader.hpp"

#include "SP3EphemerisStore.hpp"

#include "DataStructures.hpp"

#include "Counter.hpp"

#include "Epoch.hpp"


using namespace std;
using namespace gpstk;
using namespace gpstk::StringUtils;


int main(int argc, char* argv[])
{

    // conf file
    string confFileName("debug.conf");

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


    double clock_start( Counter::now() );


    SatIDSet allSatSet;

    string satListFileName("/home/kfkuang/new/ROCKET/tables/satellite_use.txt");

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
                else if(sys == "E")
                {
                    sat = SatID(id, SatID::systemGalileo);
                }
                else if(sys == "C")
                {
                    sat = SatID(id, SatID::systemBDS);
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


    CommonTime gps;

    int i(0);
    double dt(30.0);

    cout << fixed;

    // process epoch by epoch
    while( true )
    {
        gps = gps0 + i*dt;

        cout << setw(4) << CivilTime(gps).year
             << setw(4) << CivilTime(gps).month
             << setw(4) << CivilTime(gps).day
             << setw(4) << CivilTime(gps).hour
             << setw(4) << CivilTime(gps).minute
             << setprecision(3)
             << setw(10) << CivilTime(gps).second;

        for(SatIDSet::iterator it = allSatSet.begin();
            it != allSatSet.end();
            ++it)
        {
            double clkbias(0.0);

            try
            {
                clkbias = sp3Store.getXvt(*it,gps).clkbias;

                cout << setw( 5) << (*it).id
                     << setw(15) << clkbias*1e9;
            }
            catch(...)
            {
                cout << setw( 5) << (*it).id
                     << setw(15) << clkbias*1e9;
            }
        }
        cout << endl;

        i++;

        if( i*dt >= 7*24*3600.0 ) break;
    }

    double clock_end( Counter::now() );

    cerr << "time elapsed: " << setw(10) << clock_end - clock_start << endl;

    return 0;

}
