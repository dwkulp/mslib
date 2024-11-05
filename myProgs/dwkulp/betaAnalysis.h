#include <string>
#include <vector>
using namespace std;



// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input pdb
		required.push_back("pdb");

		defaultArgs.push_back("configfile");

	}

	// Storage for the vales of each optional
	string pdb;
	bool debug;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;
        vector<string> defaultArgs;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);
bool checkBetaSheet(System &_sys, int _pos);
bool neighboringStrand(System &_sys, pair<int,int> &_s1, pair<int,int> &_s2);
