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

	cout << "READ LIST"<<endl;
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

		if (line == "" || line.substr(0,1) == "#") {
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
				chainIds.push_back(line.substr(4,1));
			} else {
				cerr << "Cannot find PDB file: " << filename << endl;
			}
		} else {
			pdbs.push_back(line);
		}
	}

		
	

	fs.close();

	ChiStatistics chi;
    std::string baseMSLDir = MslTools::getBaseMSLDirectory();
    std::string chiFile = baseMSLDir + "toppar/pdb_2.3_DegOfFreedoms.txt";
    std::cout << "Chi file path: " << chiFile << std::endl;


	chi.read(chiFile);
	Transforms trans;
	AtomContainer refAtoms;
	int count = 0;
	map<string,bool> appendFile;
	for (uint i = 0; i < pdbs.size();i++){
		cout << "Working with file " << pdbs[i] ;
		if (!chainIds.empty()){
			cout << " with chainId " << chainIds[i];
		}
		cout << endl;

		// Read in PDB with no Hydrogens..
		System sys;
		if (!sys.readStructureFile(pdbs[i])) {
			cerr << "Cannot read PDB file: " << pdbs[i] << endl;
			continue;
		}

		for (uint p1 = 0; p1 < sys.positionSize(); p1++) {
			Position &pos1 = sys.getPosition(p1);
			AtomPointerVector &pos1ats = pos1.getAtomPointers();

			// Check if chainIds are not zero size and ensure the position has the same chainId as chainIds[i]
			if (!chainIds.empty() && pos1.getChainId() != chainIds[i]) {
				continue;
			}
			
			bool matchedResidueType = false;
			for (uint t = 0; t < opt.residueType.size();t++){
			  //cout << "\ttesting: "<<pos1.getResidueName()<<" "<<opt.residueType[t]<<endl;
			  if (pos1.getResidueName() == opt.residueType[t]) matchedResidueType = true;
			}
			if (!matchedResidueType) continue;

			AtomSelection sel(pos1ats);
			AtomPointerVector pos1_subSet = sel.select(opt.alignAtoms);

			vector<double> chi_angles1 = chi.getChis(pos1.getCurrentIdentity());
			cout << "CHI-RES, "<<pdbs[i]<<","<<pos1.getPositionId()<<","<< chi.getNumberChis(pos1.getCurrentIdentity());
			for (uint c = 0; c < chi.getNumberChis(pos1.getCurrentIdentity());c++){
			  cout <<","<<chi_angles1[c]<<",";
			}
			cout <<endl;
			
			if (count == 0){
				refAtoms.addAtoms(pos1_subSet);
			} else {
				if (refAtoms.size() != pos1_subSet.size()) continue;
			}
			count++;
			//cout << "FOUND A RESIDUE: "<<pdbs[i]<<" "<<pos1.toString()<<endl;

			pos1ats.saveCoor("pre");


			if (sys.positionSize() == 1) continue;

			AtomContainer neighborAtoms;
			map<string,bool> positionFound;
			for (uint a1 = 0; a1 < pos1_subSet.size();a1++){

			        //vector<int> neighbors = pos1.getCurrentIdentity().findNeighbors(opt.neighbor_dist,pos1_subSet(a1).getName(),"-N,CA,C,O,CB");
			         vector<int> neighbors = pos1.getCurrentIdentity().findNeighbors(opt.neighbor_dist,pos1_subSet(a1).getName());

				for (uint n = 0; n < neighbors.size();n++){
					Position &neighborPos = sys.getPosition(neighbors[n]);
					if (positionFound.find(neighborPos.getPositionId()) != positionFound.end()) continue;
					positionFound[neighborPos.getPositionId()] = true;
					cout << "NeighborPos: "<<neighborPos.getPositionId()<<endl;
					if (abs(neighborPos.getResidueNumber() -pos1.getResidueNumber()) < 4) continue;

					// Print Chi Angles for ASN-NAG
					if (opt.residueType[0] == "ASN" && neighborPos.getResidueName() == "NAG"){
					  vector<double> chi_angles = chi.getChis(pos1.getCurrentIdentity());
					  cout << "CHI-NAG, "<<pdbs[i]<<","<<pos1.getPositionId()<<","<<chi_angles[0]<<","<<chi_angles[1]<<endl;
					}

					stringstream ss;
					ss << opt.residueType[0]<<"_"<<neighborPos.getResidueName()<<".pdb";
					PDBWriter pout(ss.str());

					//cout << "Writing ss.str(): "<<ss.str()<<endl;
					if (appendFile.find(ss.str()) != appendFile.end()){
						pout.setOpenMode(2); // 2 == append
					}
					appendFile[ss.str()] = true;

					bool goodAlignment1 = trans.rmsdAlignment(pos1_subSet,refAtoms.getAtomPointers(),neighborPos.getAtomPointers());
					bool goodAlignment2 = trans.rmsdAlignment(pos1_subSet,refAtoms.getAtomPointers(),pos1ats);
					
					//trans.applyHistory(neighborPos.getAtomPointers());
					double rmsd = pos1_subSet.rmsd(refAtoms.getAtomPointers());
						
					if (rmsd > 1.0) {
					  pos1ats.applySavedCoor("pre");
					  continue;
					}
					//cout << "RMSDs: "<<rmsd<<" "<<" "<<goodAlignment2<<endl;

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
						
					pos1ats.applySavedCoor("pre");					

				}
			}
		}
	}
}


Options setupOptions(int theArgc, char * theArgv[]){
	Options opt;

	OptionParser OP;


	OP.setRequired(opt.required);
	OP.setAllowed(opt.optional);
	OP.readArgv(theArgc, theArgv);

	if (OP.countOptions() == 0){
		cout << "Usage:" << endl;
		cout << endl;
		cout << "superRotamerExtraction --list LIST --residueType ALA,ARG,ASN,ASP,CYS,GLN,GLU,GLY,HIS,ILE,LEU,LYS,MET,PHE,PRO,SER,THR,TRP,TYR,VAL" << endl;
		exit(0);
	}

	opt.list = OP.getString("list");
	if (OP.fail()){
		cerr << "ERROR 1111 list not specified.\n";
		exit(1111);
	}

	opt.residueType = OP.getStringVectorJoinAll("residueType");
	if (OP.fail()){
		cerr << "ERROR 1111 residueType not specified. HETERO means all non-natural amino acids\n";
		exit(1111);
	}

	opt.debug = OP.getBool("debug");
	if (OP.fail()){
	  opt.debug = false;
	}

	opt.alignAtoms = OP.getString("alignAtoms");
	if (OP.fail()){
	  opt.alignAtoms = "name N+CA+C"; // backbone atoms by default?
	}

	opt.includeAllNeighbors = OP.getBool("includeAllNeighbors");
	if (OP.fail()){
	  opt.includeAllNeighbors = false;
	}
	opt.neighbor_dist = OP.getDouble("neighborDist");
	if (OP.fail()){
	  opt.neighbor_dist = 3.5;
	}
	opt.pdbPath = OP.getString("pdbPath");
	if (OP.fail()){
		opt.pdbPath = "";
		exit(1111);
	}
	return opt;
}
