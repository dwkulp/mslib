#include <iostream>
#include <cstdlib>

#include "AtomContainer.h"
#include "AtomPointerVector.h"
#include "MslTools.h"
#include "PDBTopology.h"
#include "System.h"
#include "MslOut.h"
#include "AtomSelection.h"
#include "OptionParser.h"
#include "Transforms.h"
#include "Timer.h"
#include "ChiStatistics.h"
#include "superRotamerExtraction.h"

using namespace std;
using namespace MSL;

// MslOut 
static MslOut MSLOUT("superRotamerExtraction");

int main(int argc, char *argv[]) {

	// Read cmdline options
	Options opt = setupOptions(argc,argv);

	cout << "READ LIST: "<<opt.list<<endl;
	vector<string> pdbs;  
	vector<string> chainIds;
	ifstream fs;

	fs.open(opt.list.c_str());
	if (fs.fail()){
		cerr<<"Cannot open file "<<opt.list<<endl;
		exit(1);
	}

	bool isFirstValidLine = true;
	bool isPDBFile = false;
	while (true) {
		string line;
		getline(fs, line);

		if (fs.fail()) {
			// no more lines to read, quit the while.
			break;
		}

		if (line == "" || line[0] == '#'){
			continue;
		}
	
		// Check if the line is a PDB file or a cull-PDB type file (PDBidCHAINid) for the first non-blank line
		if (isFirstValidLine) {
			if (line.size() != 8 || line.substr(line.size() - 4) != ".pdb") {
				isPDBFile = false;
				cout << "Assuming the file type is PDB chain format (PDBidCHAINid)." << endl;
			} else {
				isPDBFile = true;
				cout << "Assuming the file type is PDB file." << endl;
			}
			isFirstValidLine = false;
		}

		if (!isPDBFile) {
			if (!MslTools::directoryExists(opt.pdbPath)) {
				cerr << "PDB path does not exist: " << opt.pdbPath << endl;
				exit(1);
			}
			string prefix = line.substr(0, 4);
			transform(prefix.begin(), prefix.end(), prefix.begin(), ::tolower);
			string filename = opt.pdbPath + "/" + prefix.substr(1, 2) + "/" + prefix + ".pdb";
			if (MslTools::fileExists(filename)) {
				pdbs.push_back(filename);
				// get chainid from line character 5 until a whitespace
				string chainId = line.substr(4, line.find(" ", 4) - 4);
				//chainIds.push_back(line.substr(4,1));
				chainIds.push_back(chainId);
			} else {
				// Get from pdb url 
				// string pdbUrl = "http://www.rcsb.org/pdb/files/" + prefix + ".pdb";
				
				stringstream ss;
				ss << MslTools::toUpper(prefix) << ".cif";
				string filename2 = ss.str();

				string pdbUrl = "https://files.rcsb.org/download/" + filename2;
				
				string command = "wget -O " + filename2 + " " + pdbUrl;
				int retVal = system(command.c_str());
				
				if (MslTools::fileExists(filename2)){
					pdbs.push_back(filename2);
					string chainId = line.substr(4, line.find(" ", 4) - 4);
					chainIds.push_back(chainId);
					
					//chainIds.push_back(line.substr(4,1));
				} else {
					cerr << "Cannot find get PDB file: " << filename << " or remote get file: "<<filename2<<endl;
				}
			}
		} else {
			pdbs.push_back(line);
		}
	}

		
	

	fs.close();


	ChiStatistics chi;

	
	if (opt.chi){
		std::string baseMSLDir = MslTools::getBaseMSLDirectory();
    	std::string chiFile = baseMSLDir + "toppar/pdb_2.3_DegOfFreedoms.txt";
	    std::cout << "Chi file path: " << chiFile << std::endl;
		if (!MslTools::fileExists(chiFile)){
			cerr << "Cannot find chi file: "<<chiFile<<endl;
			exit(1);
		}
		chi.read(chiFile);
	}
	Transforms trans;
	AtomContainer refAtoms;
	int count = 0;
	map<string,bool> appendFile;
	int residueTypeCount = 0;
	int residueTypeCountWithNeighbors = 0;

	for (uint i = 0; i < pdbs.size();i++){
		cout << "Working with file " << pdbs[i] ;
		if (!chainIds.empty()){
			cout << " with chainId " << chainIds[i] << ".";
		}
		// Add endl if in debug because more things will be printed out
		if (opt.debug)
			cout << endl;

		// Read in PDB with no Hydrogens..
		System sys;
		if (!sys.readStructureFile(pdbs[i])) {
			cerr << "Cannot read PDB file: " << pdbs[i] << endl;
			continue;
		}

		// Only count each pair once
		map<string,bool> positionPairFound;
		int residueTypeCountWithNeighborsPerPDB = 0;
		sys.getAtomPointers().saveCoor("pre");

		for (uint p1 = 0; p1 < sys.positionSize(); p1++) {
			Position &pos1 = sys.getPosition(p1);
			AtomPointerVector &pos1ats = pos1.getAtomPointers();

			if (opt.debug > 1 ){
				cout << "Position: "<<pos1.getPositionId()<<endl;
			}

			// Check if chainIds are not zero size and ensure the position has the same chainId as chainIds[i]
			if (!chainIds.empty() && pos1.getChainId() != chainIds[i]) {
				continue;
			}

			if (opt.debug){
				cout << "Chain Found Position: "<<pos1.getPositionId()<<endl;
			}
			bool matchedResidueType = false;
			for (uint t = 0; t < opt.residueType.size();t++){
			  //cout << "\ttesting: "<<pos1.getResidueName()<<" "<<opt.residueType[t]<<endl;
			  if (pos1.getResidueName() == opt.residueType[t]) {
				matchedResidueType = true;
            	residueTypeCount++;
			  }
        	}


			if (!matchedResidueType) continue;
			if (opt.debug){
				cout << "Residue Type Found Position: "<<pos1.getPositionId()<<endl;
			}
			AtomSelection sel(pos1ats);
			AtomPointerVector pos1_subSet;
			AtomPointerVector pos1_neighborTestingSet;
			if (opt.alignAtoms != ""){
				pos1_subSet = sel.select(opt.alignAtoms);
				if (opt.debug){
					cout << "Align Atoms: "<<pos1_subSet.size()<<endl;
				}
			}
			if (opt.neighborTestAtoms != ""){
				pos1_neighborTestingSet = sel.select(opt.neighborTestAtoms);
				if (opt.debug){
					cout << "Neighbor Testing Atoms: "<<pos1_neighborTestingSet.size()<<endl;
				}
			} else {
				if (opt.alignAtoms != ""){
					pos1_neighborTestingSet = sel.select(opt.alignAtoms);
				} 
			}


			if (opt.chi){
				vector<double> chi_angles1 = chi.getChis(pos1.getCurrentIdentity());
				cout << "CHI-RES, "<<pdbs[i]<<","<<pos1.getPositionId()<<","<< chi.getNumberChis(pos1.getCurrentIdentity());
				for (uint c = 0; c < chi.getNumberChis(pos1.getCurrentIdentity());c++){
				  cout <<","<<chi_angles1[c]<<",";
				}
				cout <<endl;
			}
			
			if (count == 0){
				refAtoms.addAtoms(pos1_subSet);
			} else {
				if (refAtoms.size() != pos1_subSet.size()) continue;
			}
			count++;
			//cout << "FOUND A RESIDUE: "<<pdbs[i]<<" "<<pos1.toString()<<endl;

			pos1ats.saveCoor("pre");


			if (sys.positionSize() == 1) continue;

			if (opt.debug){
				cout << "Valid Position: "<<pos1.getPositionId()<<endl;
			}

			AtomContainer neighborAtoms;
			map<string,bool> positionFound;
			for (uint a1 = 0; a1 < pos1_neighborTestingSet.size();a1++){

			        //vector<int> neighbors = pos1.getCurrentIdentity().findNeighbors(opt.neighbor_dist,pos1_subSet(a1).getName(),"-N,CA,C,O,CB");
			    
				// Select all atoms using this pos1_neighborTestingSet atom to look for neighbors
				// if no neighborTestAtoms are specified, then use all atoms
				vector<int> neighbors;
				if (opt.neighborTestAtoms == ""){
					neighbors = pos1.getCurrentIdentity().findNeighbors(opt.neighbor_dist,pos1_neighborTestingSet(a1).getName());
				} else {
					neighbors = pos1.getCurrentIdentity().findNeighbors(opt.neighbor_dist,pos1_neighborTestingSet(a1).getName(), pos1_neighborTestingSet(a1).getName());
				}
				
				if (opt.debug > 1){
					cout << "Neighbors: "<<neighbors.size()<<endl;
				}

				for (uint n = 0; n < neighbors.size();n++){
					Position &neighborPos = sys.getPosition(neighbors[n]);
					if (positionFound.find(neighborPos.getPositionId()) != positionFound.end()) continue;
					positionFound[neighborPos.getPositionId()] = true;
					

					if (positionPairFound.find(pos1.getPositionId() + neighborPos.getPositionId()) != positionPairFound.end() || positionPairFound.find(neighborPos.getPositionId()+ pos1.getPositionId()) != positionPairFound.end()) continue;

					positionPairFound[pos1.getPositionId() + neighborPos.getPositionId()] = true;
					positionPairFound[neighborPos.getPositionId() + pos1.getPositionId()] = true;

					if (opt.debug) {
						cout << "Position1: " << pos1.getPositionId() <<" "<<pos1.getResidueName()<<" NeighborPos: "<<neighborPos.getPositionId()<<" "<<neighborPos.getResidueName()<<endl;
					}

					// Don't get local neighbors
					if (abs(neighborPos.getResidueNumber() - pos1.getResidueNumber()) < 4) continue;

					// Check if neighborPos.getResidueName() is in partnerResidueType
		            if (!opt.partnerResidueType.empty() && std::find(opt.partnerResidueType.begin(), opt.partnerResidueType.end(), neighborPos.getResidueName()) == opt.partnerResidueType.end()) {
        		        continue;
            		}

					if (opt.debug){
						cout << "Position2: " << pos1.getPositionId() <<" "<<pos1.getResidueName()<<" NeighborPos: "<<neighborPos.getPositionId()<<" "<<neighborPos.getResidueName()<<endl;
					}

					// Print Chi Angles for ASN-NAG
					if (opt.residueType[0] == "ASN" && neighborPos.getResidueName() == "NAG"){
					  vector<double> chi_angles = chi.getChis(pos1.getCurrentIdentity());
					  cout << "CHI-NAG, "<<pdbs[i]<<","<<pos1.getPositionId()<<","<<chi_angles[0]<<","<<chi_angles[1]<<endl;
					}
					residueTypeCountWithNeighbors++;
					residueTypeCountWithNeighborsPerPDB++;	

					stringstream ss;
					ss << MslTools::getFileName(opt.list) << "_" << opt.residueType[0]<<"_"<<neighborPos.getResidueName()<<".pdb";
					PDBWriter pout(ss.str());

					if (opt.debug){
						cout << "Writing ss.str(): "<<ss.str()<<endl;
					}

					if (appendFile.find(ss.str()) != appendFile.end()){
						pout.setOpenMode(2); // 2 == append
					}
					appendFile[ss.str()] = true;

					if (opt.alignAtoms != ""){
						bool goodAlignment1 = trans.rmsdAlignment(pos1_subSet,refAtoms.getAtomPointers(),neighborPos.getAtomPointers());
						bool goodAlignment2 = trans.rmsdAlignment(pos1_subSet,refAtoms.getAtomPointers(),pos1ats);
					
						//trans.applyHistory(neighborPos.getAtomPointers());
						double rmsd = pos1_subSet.rmsd(refAtoms.getAtomPointers());
						
						if (rmsd > 1.0) {
					  		pos1ats.applySavedCoor("pre");
					 		continue;
						}
						if (opt.debug) {
							cout << "RMSDs: "<<rmsd<<" "<<" "<<goodAlignment2<<endl;
						}
					}

					AtomPointerVector ats;
					ats += pos1ats;
					ats += neighborPos.getAtomPointers();
					
					if (opt.includeAllNeighbors){
					  for (uint m = 0; m < neighbors.size();m++){
					    if (m != n){
					      ats += sys.getPosition(neighbors[m]).getAtomPointers();
					    }
					  }
					}
					pout.open();
					pout.write(ats,false,false,true);
					pout.close();

					// Reset the coordinates
					if (opt.alignAtoms != ""){
						sys.getAtomPointers().applySavedCoor("pre");
					}

				}
			}
		}
		cout << " pairs of residues found: "<<residueTypeCountWithNeighborsPerPDB<<endl;
	}
	cout << "Total number of residues that are of specified type: "<< residueTypeCount << " with "<< residueTypeCountWithNeighbors<<  " neighbors in " << pdbs.size()<< " pdb files."<<endl;
}


Options setupOptions(int argc, char* argv[]) {
    Options opt;

    OptionParser OP;
    OP.setRequired(opt.required);
    OP.setAllowed(opt.optional);
    //OP.setDefaultArguments(opt.defaultArgs); // a pdb file value can be given as a default argument without the --pdbfile option
    OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"
    OP.readArgv(argc, argv);

    if (OP.countOptions() == 0) {
        std::cout << "Usage:" << std::endl;
        std::cout << "  superRotamerExtraction --pdbPath <PDB_PATH> [options]" << std::endl;
        std::cout << std::endl;
        std::cout << "Required options:" << std::endl;
        for (uint i = 0; i < opt.required.size(); i++) {
            std::cout << "  --" << opt.required[i] << std::endl;
        }
        std::cout << std::endl;
        std::cout << "Optional options:" << std::endl;
        for (uint i = 0; i < opt.optional.size(); i++) {
            std::cout << "  --" << opt.optional[i] << std::endl;
        }
        std::cout << std::endl;
        std::exit(0);
    }
	opt.list = OP.getString("list");
	if (OP.fail()){
		std::cerr << "ERROR: list not specified.\n";
		std::exit(1);
	}
    opt.pdbPath = OP.getString("pdbPath");
    if (OP.fail()) {
		opt.pdbPath = "";
    }

    opt.residueType = OP.getStringVectorJoinAll("residueType");
	if (OP.fail()) {
        std::cerr << "WARNING: residueType not specified.\n";
	}

    opt.partnerResidueType = OP.getStringVectorJoinAll("partnerResidueType");
    if (OP.fail()) {
        opt.partnerResidueType.clear(); // Default to empty vector if not specified
    }

	// Handle other options...
	opt.alignAtoms = OP.getString("alignAtoms");
	if (OP.fail()) {
		opt.alignAtoms = "name N+CA+C"; // backbone atoms by default?
	}

	opt.neighborTestAtoms = OP.getString("neighborTestAtoms");
	if (OP.fail()) {
		opt.neighborTestAtoms = "";
	}

	opt.includeAllNeighbors = OP.getBool("includeAllNeighbors");
	if (OP.fail()) {
		opt.includeAllNeighbors = false;
	}
	opt.neighbor_dist = OP.getDouble("neighborDist");
	if (OP.fail()) {
		opt.neighbor_dist = 2.5;
	}
	opt.chi = OP.getBool("chi");
	if (OP.fail()) {
		opt.chi = false;
	}

    opt.debug = OP.getInt("debug");
    if (OP.fail()){
      opt.debug = 0;
	}

    return opt;
}
