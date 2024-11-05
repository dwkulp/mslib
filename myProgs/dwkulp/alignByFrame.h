#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input residues
		required.push_back("pdb1");
		required.push_back("pdb2");
		
		// Selection
		required.push_back("sel1");
		required.push_back("sel2");


	}

	// Storage for the vales of each optional
        string pdb1;
        string pdb2;
	string sel1;	
	string sel2;	

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

