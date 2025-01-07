/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)n
 Copyright (C) 2009 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan

This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, 
 USA, or go to http://www.gnu.org/copyleft/lesser.txt.
----------------------------------------------------------------------------
*/
#include <string>
#include <map>
#include <fstream>
#include <signal.h>
#include "OptionParser.h"
#include "Timer.h"
#include "System.h"
#include "AtomSelection.h"
#include "SasaCalculator.h"
#include "AtomSelection.h"
#include "getDeltaDistSASA.h"
#include "Transforms.h"

using namespace std;

using namespace MSL;

// Global objects
Timer t;
double startTime = 0.0;

/**
 * @file getDeltaDistSASA.cpp
 * @brief This program calculates the difference in Solvent Accessible Surface Area (SASA) between two protein structures.
 *
 * The program reads two PDB files, calculates the SASA for each residue in both structures, and then computes the difference
 * in SASA between the corresponding residues of the two structures. It also normalizes the SASA values based on reference
 * SASA values for each amino acid type.
 *
 * The program outputs the position ID, distance between C-alpha atoms, SASA values for both structures, delta SASA, 
 * normalized SASA values, and delta normalized SASA for each residue.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 * @return int Exit status of the program.
 */
int main(int argc, char *argv[]){

	// Option Parser
	Options opt = setupOptions(argc, argv);

	// Parse epitopes and map residue numbers to epitope names
	map<string, string> epitopeMap;
	map<string, vector<int> > epitopeMapBack;
	for (const string &epitope : opt.epitopes) {
		size_t colonPos = epitope.find(':');
		if (colonPos != string::npos) {
			string name = epitope.substr(0, colonPos);
			string positions = epitope.substr(colonPos + 1);
			size_t dashPos = positions.find('-');
			if (dashPos != string::npos) {
				string resNum1 = positions.substr(0, dashPos);
				string resNum2 = positions.substr(dashPos + 1);
				int start = stoi(resNum1);
				int end = stoi(resNum2);
				for (int resNum = start; resNum <= end; ++resNum) {
					epitopeMap[to_string(resNum)] = name;
					epitopeMapBack[name].push_back(resNum);
				}
			}
		}
	}

	// Create a system from the structural input options
	startTime = t.getWallTime();

	// Read in PDB1
	System sys1;
	sys1.readStructureFile(opt.pdb1);

	// Read in PDB2
	System sys2;
	sys2.readStructureFile(opt.pdb2);

	// Compute geometric center for posId defined by referencePoint for each sys1 and sys2
	CartesianPoint refPoint1, refPoint2;
	if (opt.referencePoint != "") {
		AtomSelection sel1(sys1.getAtomPointers());
		AtomSelection sel2(sys2.getAtomPointers());

		AtomPointerVector atoms1 = sel1.select(opt.referencePoint);
		AtomPointerVector atoms2 = sel2.select(opt.referencePoint);

		if (atoms1.size() == 0 || atoms2.size() == 0 || atoms1.size() != atoms2.size()) {
			cerr << "ERROR: referencePoint not found in one of the structures.\n";
			exit(1111);
		}

		refPoint1 = atoms1.getGeometricCenter();
		refPoint2 = atoms2.getGeometricCenter();

		// Align sys2 to sys1 by refPoint2 to refPoint1
		Transforms tm;
		if (!tm.rmsdAlignment(atoms2, atoms1, sys2.getAtomPointers())) {
			cerr << "Alignment failed!" << endl;
			exit(1);
		}

		//double rmsd = tm.getRMSD();
		double rmsd = atoms1.rmsd(atoms2);
		printf("RMSD: %f\n", rmsd);
	}


     /*
	SASA reference:
	Protein Engineering vol.15 no.8 pp.659–667, 2002
	Quantifying the accessible surface area of protein residues in their local environment
	Uttamkumar Samanta Ranjit P.Bahadur and  Pinak Chakrabarti
      */
 	 map<string,double> refSasa;
	refSasa["G"] = 83.91;
	refSasa["A"] = 116.40;
	refSasa["S"] = 125.68;
	refSasa["C"] = 141.48;
	refSasa["P"] = 144.80;
	refSasa["T"] = 148.06;
	refSasa["D"] = 155.37;
	refSasa["V"] = 162.24;
	refSasa["N"] = 168.87;
	refSasa["E"] = 187.16;
	refSasa["Q"] = 189.17;
	refSasa["I"] = 189.95;
	refSasa["L"] = 197.99;
	refSasa["H"] = 198.51;
	refSasa["K"] = 207.49;
	refSasa["M"] = 210.55;
	refSasa["F"] = 223.29;
	refSasa["Y"] = 238.30;
	refSasa["R"] = 249.26;
	refSasa["W"] = 265.42;

	// Identify Exposed positions
	SasaCalculator scalc1(sys1.getAtomPointers());
	scalc1.calcSasa();

	SasaCalculator scalc2(sys2.getAtomPointers());
	scalc2.calcSasa();
  
	// Open outputfile
	ofstream outFile;
	outFile.open("dist_sasa.csv");

	// print formated header to outfile
  	outFile << MslTools::stringf("%-10s %-4s %10s %10s %10s %10s %10s %10s %10s %10s %s\n",
					   "PosID", "Res", "Dist", "DeltaDist", "SASA1", "SASA2", "DeltaSASA", "NormSASA1", "NormSASA2", "DeltaNormSASA", "Epitope");
	for (uint i = 0; i < sys1.positionSize(); i++) {
		string posId = sys1.getPosition(i).getPositionId();

		if (sys2.positionExists(posId) && sys1.getPosition(posId).getResidueName() == sys2.getPosition(posId).getResidueName()) {
			Position &pos1 = sys1.getPosition(posId);
			Position &pos2 = sys2.getPosition(posId);

			double dist = 0.0;
			if (pos1.atomExists("CA") && pos2.atomExists("CA")) {
				dist = pos1.getAtom("CA").distance(pos2.getAtom("CA"));
			}

			double deltaDeltaDist = 0.0;
			if (opt.referencePoint != "") {
				double deltaDist1 = 0.0;
				double deltaDist2 = 0.0;
				if (pos1.atomExists("CA")) {
					deltaDist1 = pos1.getAtom("CA").getCoor().distance(refPoint1);
				}
				if (pos2.atomExists("CA")) {
					deltaDist2 = pos2.getAtom("CA").getCoor().distance(refPoint2);
				}
				deltaDeltaDist = deltaDist1 - deltaDist2;
			}
			
			double sasa1 = scalc1.getResidueSasa(posId);
			double sasa2 = scalc2.getResidueSasa(posId);

			string resAA1 = MslTools::getOneLetterCode(pos1.getResidueName());
			string resAA2 = MslTools::getOneLetterCode(pos2.getResidueName());

			double normSasa1 = sasa1 / refSasa[resAA1];
			double normSasa2 = sasa2 / refSasa[resAA2];
			double deltaNormSasa = normSasa1 - normSasa2;
			double deltaSasa = sasa1 - sasa2;

			string resNumStr = to_string(pos1.getResidueNumber());
			string epitopeName = epitopeMap.count(resNumStr) ? epitopeMap[resNumStr] : "None";

	

			outFile << MslTools::stringf("%-10s %-4s %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %s\n",
				   posId.c_str(), pos1.getResidueName().c_str(), dist, deltaDeltaDist, sasa1, sasa2, deltaSasa, normSasa1, normSasa2, deltaNormSasa, epitopeName.c_str());
		} else {
			cerr << "ERROR: Position " << posId << " not found in second structure.\n";
		}
	}
	outFile.close();

	// Figure out which epitopes are interacting in sys1 or sys2 (note that chain 'A' is hardcoded here)
	sys1.getChains()[0]->setChainId("A");
	sys2.getChains()[0]->setChainId("A");
	for (const auto &epitope1 : epitopeMapBack) {
		for (const auto &epitope2 : epitopeMapBack){
			if (epitope1.first != epitope2.first){	
				//cout << "Epitope1: " << epitope1.first << ", Epitope2: " << epitope2.first << endl;
				// measure all pairwise distances between epitope1 and epitope2 residues in sys1 and sys2
				int countSys1 = 0;
				int countSys2 = 0;
				
				for (int resNum1 : epitope1.second) {
					for (int resNum2 : epitope2.second) {
						
						bool residuesInContact1 = false;
						bool residuesInContact2 = false;

						if (sys1.positionExists("A",resNum1,"") && sys1.positionExists("A",resNum2,"")) {
							Position &pos1 = sys1.getPosition("A",resNum1,"");
							Position &pos2 = sys1.getPosition("A", resNum2,"");
							if (pos1.atomExists("CA") && pos2.atomExists("CA")) {
								double dist = pos1.getAtom("CA").distance(pos2.getAtom("CA"));
								if (dist < 8.0) {
									countSys1++;
									residuesInContact1 = true;
								}
								
							}
						}
						if (sys2.positionExists("A",resNum1,"") && sys2.positionExists("A",resNum2,"")) {
							Position &pos1 = sys2.getPosition("A",resNum1,"");
							Position &pos2 = sys2.getPosition("A",resNum2,"");
							if (pos1.atomExists("CA") && pos2.atomExists("CA")) {
								double dist = pos1.getAtom("CA").distance(pos2.getAtom("CA"));
								if (dist < 8.0) {
									countSys2++;
									residuesInContact2 = true;
								}

							}
						}

						if ((residuesInContact1 && !residuesInContact2) || (!residuesInContact1 && residuesInContact2)){
							printf("Epitope pair (%-5s, %-5s): Residues %-5d and %-5d change contact\n", 
								   epitope1.first.c_str(), epitope2.first.c_str(), resNum1, resNum2);
						}
					}
				}
				if ( abs(countSys1 - countSys2) > 9){
					printf("Epitope pair (%-5s, %-5s): Sys1 count = %-5d, Sys2 count = %-5d\n", 
						   epitope1.first.c_str(), epitope2.first.c_str(), countSys1, countSys2);
				}

			}

		}
	}
	// Figure out which epitopes are interacting in sys1 or sys2 (note that chain 'A' is hardcoded here) keep track of specific residue-residue contacts that change between sys1 and sys2 store a hash of resnum-resnum and epitope
	

	
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
		cout << "getDeltaDistSASA --pdb1 PDB1 --pdb2 PDB2 --epitopes epitope1:resNum1-resNum2,epitope2:resNum3-resNum4 --referencePoint referencePoint" << endl;
		exit(0);
	}
	opt.pdb1 = OP.getString("pdb1");
	if (OP.fail()){
		cerr << "ERROR: pdb1 not specified.\n";
		exit(1111);
	}
	opt.pdb2 = OP.getString("pdb2");
	if (OP.fail()){
		cerr << "ERROR: pdb2 not specified.\n";
		exit(1111);
	}
	opt.epitopes = OP.getStringVectorJoinAll("epitopes");
	if (OP.fail()){
		cerr << "ERROR: epitopes not specified.\n";
		exit(1111);
	}

	opt.referencePoint = OP.getString("referencePoint");
	if (OP.fail()){
		cerr << "ERROR: referencePoint not specified.\n";
		exit(1111);
	}
	
	return opt;
}







