#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// List of PDBs to search against
		required.push_back("pdb");
		required.push_back("atom1");
		required.push_back("atom2");
		required.push_back("atom1sel");
		required.push_back("atom2sel");
		required.push_back("terminiDistance");
		
		optional.push_back("debug");

	}

	// Storage for the vales of each optional
	string pdb;
        string atom1;
        string atom2;
        string atom1sel;
        string atom2sel;
        bool terminiDist;
	bool debug;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

