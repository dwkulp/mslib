#include <string>
#include <vector>
#include <map>

using namespace std;



// Define Input Options and Store them.
struct Options {

  // Set up options here...
  Options(){

    // Input pdb
    required.push_back("pdb");
    optional.push_back("sasa");
    optional.push_back("pick_random_SSEs");
    optional.push_back("include_loops");
    optional.push_back("neighbor_dist");
    optional.push_back("num_contacting");
    optional.push_back("msa"); // Added msa option
    optional.push_back("aa_freq_threshold"); // Added aa_freq option
    optional.push_back("exposed_pos_threshold"); // Added exposed_pos_threshold option
    
  }

  // Storage for the values of each optional
  string pdb;
  bool sasa;
  bool include_loops;
  int pick_random_SSEs;
  double neighbor_dist;
  int num_contacting;
  string msa; // Added msa storage
  double aa_freq_threshold; // Added aa_freq storage
  int exposed_pos_threshold; // Added exposed_pos_threshold storage
  
  // Storage for different types of options
  vector<string> required;
  vector<string> optional;
  vector<string> defaultArgs;

};

struct sasaData {
  int total_pos;
  int exposed_pos;
  double total_nSasa;

  sasaData(){
    total_pos = 0;
    exposed_pos = 0;
    total_nSasa = 0.0;
  }

   sasaData(const sasaData &_rhs){
     total_pos      = _rhs.total_pos;
     exposed_pos    = _rhs.exposed_pos;;
     total_nSasa    = _rhs.total_nSasa;
  }
  void operator=(const sasaData &_rhs){
    total_pos      = _rhs.total_pos;
     exposed_pos    = _rhs.exposed_pos;;
     total_nSasa    = _rhs.total_nSasa;
  }
  
};

/*
	SASA reference:
	Protein Engineering vol.15 no.8 pp.659–667, 2002
	Quantifying the accessible surface area of protein residues in their local environment
	Uttamkumar Samanta Ranjit P.Bahadur and  Pinak Chakrabarti
*/
map<string,double> refSasa;



// Helper functions to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);
bool checkBetaSheet(System &_sys, int _pos);
vector<map<string, int> > computeAADiversity(System &sys, const string &msaFilePath, const string &pdbFileName);




