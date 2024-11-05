#include <vector>
using namespace MSL;
using namespace std;
struct Options {

	// Set up options here...
	Options(){

		// PDB list
		required.push_back("pdb");
		required.push_back("disulfPdb");
		optional.push_back("positions");
		optional.push_back("tol");
		optional.push_back("modelBest");
		optional.push_back("writeOutAll");
		optional.push_back("skip_local_disulfides");
		optional.push_back("only_inter_chain");
		optional.push_back("fasta");
		optional.push_back("appendName");
		optional.push_back("binaryDB");
		

	}




	// Storage for the vales of each option
	string pdb;
        string disulfPdb;
        double tol;
        bool modelBest;
        bool writeOutAll;
        int skip_local_disulfides;
        bool only_inter_chain;
        bool fasta;
        string appendName;
        bool binaryDB;
        vector<string> specific_positions;

	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};


Options setupOptions(int theArgc, char * theArgv[]);

