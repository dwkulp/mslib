#include <vector>
using namespace MSL;
using namespace std;
struct Options {

	// Set up options here...
	Options(){

		// PDB list
		required.push_back("pdblist");
		required.push_back("tooManyClashes");
		required.push_back("resTypes");
		required.push_back("dist");
		required.push_back("reportClashes");
	}




	// Storage for the vales of each option
	string pdblist;
        int tooManyClashes;
        vector<string> residueTypes;
        int dist;
        bool printClashes;

	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};


Options setupOptions(int theArgc, char * theArgv[]);

