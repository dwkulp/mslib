#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

    // Set up options here...
    Options(){

        // List of PDBs to search against
        required.push_back("pdb1");
        required.push_back("pdb2");
        
        optional.push_back("debug");
        optional.push_back("epitopes");
        optional.push_back("referencePoint");

    }

    // Storage for the values of each optional
    string pdb1;
    string pdb2;
    bool debug;
    vector<string> epitopes;
    string referencePoint;

    // Storage for different types of options
    vector<string> required;
    vector<string> optional;

};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

