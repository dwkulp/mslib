#include <string>
#include <vector>
using namespace std;



// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input pdb
		required.push_back("pdb");
		optional.push_back("numModels");


	}

	// Storage for the vales of each optional
	string pdb;
        int numModels;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;
        vector<string> defaultArgs;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);


