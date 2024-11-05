#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input residues
		required.push_back("pdb");
		required.push_back("sele1");
		required.push_back("sele2");
	}

	// Storage for the vales of each optional
	string pdb;
  	string sele1;
    	string sele2;

  
	// Storage for different types of options
	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);
void translate(AtomPointerVector &_ats, double _dist, CartesianPoint _center, CartesianPoint _GC, Line _line);
