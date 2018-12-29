#pragma ident "$Id$"

#include <algorithm>

#include "ConfDataReader.hpp"

#include "Rinex3ClockStream.hpp"
#include "Rinex3ClockHeader.hpp"
#include "Rinex3ClockData.hpp"

#include "DataStructures.hpp"

#include "Counter.hpp"

#include "Epoch.hpp"


using namespace std;
using namespace gpstk;
using namespace gpstk::StringUtils;


int main(int argc, char* argv[])
{
    if(argc < 1) exit(-1);

    string station( argv[1] );

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


    vector<string>::const_iterator clk_it = clkFileListVec.begin();
    while(clk_it != clkFileListVec.end())
    {
        clkFile = (*clk_it);

        Rinex3ClockStream strm(clkFile.c_str());

        strm.exceptions(std::ios::failbit);


        Rinex3ClockHeader header;
        Rinex3ClockData data;

        try
        {
            strm >> header;

            while(strm >> data)
            {
                if(data.datatype == string("AR"))
                {
                    if(data.site == station)
                    {
                        cout << fixed << setprecision(3)
                             << setw(10) << data.site
                             << setw(5) << CivilTime(data.time).year
                             << setw(5) << CivilTime(data.time).month
                             << setw(5) << CivilTime(data.time).day
                             << setw(5) << CivilTime(data.time).hour
                             << setw(5) << CivilTime(data.time).minute
                             << setw(10) << CivilTime(data.time).second
                             << setw(20) << data.bias*1e9
                             << setw(20) << data.drift*1e9
                             << endl;
                    }
                }
            }
        }
        catch(...)
        {
            ++clk_it;
            continue;
        }

        ++clk_it;

        strm.close();
    }

    return 0;
}
