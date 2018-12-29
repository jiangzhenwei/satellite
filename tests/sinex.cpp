
#include "SinexStream.hpp"
#include "SinexHeader.hpp"
#include "SinexData.hpp"

using namespace std;
using namespace gpstk;

int main(int argc, char *argv[])
{
    if(argc < 2)
    {
        cerr << "No Sinex file input." << endl;
        exit(-1);
    }

    try
    {
        Sinex::Data data;

        cout << "Reading " << argv[1] << "..." << endl;
        Sinex::Stream input(argv[1]);
        input.exceptions(fstream::eofbit | fstream::failbit);
        input >> data;
        cout << "Done." << endl;

//        data.dump(cout);

//        cout << "Writing data to sinex_test.out ..." << endl;
//        Sinex::Stream output("sinex_test.out", ios::out | ios::ate);
//        output.exceptions(fstream::eofbit | fstream::failbit);
//        output << data;
//        cout << "Done." << endl;

        exit(0);
    }
    catch(Exception& e)
    {
        cerr << e;
    }
    catch(...)
    {
        cerr << "processing error." << endl;
    }

    exit(1);

}
