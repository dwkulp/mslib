#include <vector>
#include "PyMolVisualization.h"
struct Options {

	// Set up options here...
	Options(){

		// PDB
		required.push_back("pdb");
		required.push_back("fragdb");
		optional.push_back("regex");
		optional.push_back("sel");
		optional.push_back("rmsd");
		optional.push_back("dump");
		optional.push_back("outpdb");
		optional.push_back("maxMatches");

	}

	// Storage for the vales of each option
	string pdb;
        string fragdb;
        string regex;
        string sel;
        double rmsd;
        int dump;
        string outpdb;
        string pdbDir;
        bool includeFullFile;
        int maxMatches;

	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};

Options setupOptions(int theArgc, char * theArgv[]);

