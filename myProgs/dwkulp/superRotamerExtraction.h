#include <string>
#include <vector>
using namespace std;

// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input residues
		required.push_back("list");
		required.push_back("residueType");
		required.push_back("separateByResidueType");

		// Number of rotamers
		optional.push_back("debug");
		optional.push_back("alignAtoms");
		optional.push_back("includeAllNeighbors");
		optional.push_back("neighborDist");
		optional.push_back("pdbPath");
		optional.push_back("partnerResidueType");
		optional.push_back("chi");
		optional.push_back("neighborTestAtoms");

	}

	// Storage for the values of each optional
	string list;
	vector<string> residueType;
	string alignAtoms;
	string neighborTestAtoms;
	bool separateByResidueType;
	bool includeAllNeighbors;
	bool chi;
	int debug;
	double neighbor_dist;
	string pdbPath;
	vector<string> partnerResidueType;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;

};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

