#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input residues
	        required.push_back("pdbs");
		
		// 
		optional.push_back("numChains");
		optional.push_back("numCloseAts");
		optional.push_back("distCutoff");

		// Config file
	        defaultArgs.push_back("configfile");

	}

	// Storage for the vales of each optional
	string pdbs;
        int numChains;
        int numCloseAts;
        double distCutoff;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

