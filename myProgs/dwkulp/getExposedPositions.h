#ifndef GET_EXPOSED_POSITIONS_H
#define GET_EXPOSED_POSITIONS_H

#include <vector>


struct Options {

	// Set up options here...
	Options(){

		// Input PDB File
		required.push_back("pdb");
		optional.push_back("selection");
		optional.push_back("sasaCutoff");
		optional.push_back("pymol");

	}

    string pdbFile;
    double percentSasa;
    string selectPositions;
    bool pymol;

    std::vector<std::string> required;
	std::vector<std::string> optional;
};

struct subSeq {
    string sse;
    string seq;
    vector<int> positions;
};

// Function to get exposed positions
vector<subSeq> getExposedPositions(System &_sys, Options &_opt);
Options setupOptions(int theArgc, char * theArgv[]);

#endif // GET_EXPOSED_POSITIONS_H