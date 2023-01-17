#include <vector>
using namespace MSL;
using namespace std;
struct Options {

	// Set up options here...
	Options(){

		// PDB list
		required.push_back("pdb");
		required.push_back("seq");
		required.push_back("table");
		required.push_back("out");
		required.push_back("fasta");
	}




	// Storage for the vales of each option
  string pdb;
  string seq;
  string table;
  string out;
  string fasta;

	vector<string> required;
	vector<string> optional;
	vector<string> defaultArgs;
};


Options setupOptions(int theArgc, char * theArgv[]);

