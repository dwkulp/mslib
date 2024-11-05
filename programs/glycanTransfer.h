#include <string>
#include <vector>

using namespace std;

struct DoF {
        vector<string> residueNames; // Could be 1 residue or 2 residues (for a dihedral between two residues).
	vector<string> atomNames;
	double value;
	string type; // bond, angle, dihedral, improper 
        double range_start;
        double range_stop;
        double range_step;

        DoF(){
	  range_start = 0;
	  range_stop  = 360;
	  range_step  = 5;
	}
};


// Define Input Options and Store them.
struct Options {

	// Set up options here...
	Options(){

		// Input pdb
		required.push_back("pdb");
		required.push_back("newRes");
		required.push_back("positions");

		// Input fasta
		required.push_back("dihedrals");
		optional.push_back("range");

		optional.push_back("clashTol");
		optional.push_back("alignAtoms");
		optional.push_back("rmsdAtoms");
		optional.push_back("outPdb");
		optional.push_back("randomize");
		optional.push_back("maxIter");
		optional.push_back("doffile");
		optional.push_back("configfile");
		defaultArgs.push_back("configfile");

	}

	// Storage for the vales of each optional
	string pdb;
        string newRes;
        string outPdb;
        int randomize;
        int maxIter;
        vector<string> positions; 
        vector<vector<string> > dihedrals;
        vector<string> range;
        string alignAtoms;
        string rmsdAtoms;
        int clashTol;
        string configfile;
        string doffile;
         vector<DoF> DoFs;


	// Storage for different types of options
	vector<string> required;
	vector<string> optional;
        vector<string> defaultArgs;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);
int countClashes(AtomPointerVector &_ats, AtomPointerVector &_b,bool _checkInternal=true, int _clashTol=100);
int countClashes(AtomPointerVector &_ats, Atom3DGrid &_grid, bool _checkInternal=true, int _clashTol=100);
void setDihedral(System &_sys, DoF &_dih,double _value);
