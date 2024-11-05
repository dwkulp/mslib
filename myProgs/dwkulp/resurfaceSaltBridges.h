#include <string>
#include <vector>
using namespace std;



// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input pdb
		required.push_back("pdb");
		required.push_back("topfile");
		required.push_back("parfile");
		required.push_back("rotlib");
		required.push_back("numStructuralModels");
		required.push_back("numSequenceModels");
		
		optional.push_back("numMCcycles");
		optional.push_back("startMCtemp");
		optional.push_back("endMCtemp");
		optional.push_back("sb_prop_table");
		optional.push_back("scoreOnly");
		optional.push_back("selectPositions");
		optional.push_back("percentSasa");
		optional.push_back("fasta");
		optional.push_back("appendName");
		defaultArgs.push_back("configfile");

	}

	// Storage for the vales of each optional
	string pdb;
        string sb_prop_table;
        string topfile;
        string parfile;
        string rotlib;
        int numStructuralModels;
        int numSequenceModels;
        int numMCcycles;
        double startMCtemp;
        double endMCtemp;
        string selectPositions;
        double percentSasa;
        bool scoreOnly;
        bool fasta_only;
        string appendName;
	bool debug;

	// Storage for different types of options
	vector<string> required;
	vector<string> optional;
        vector<string> defaultArgs;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);

