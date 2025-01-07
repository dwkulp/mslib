/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
 Sabareesh Subramaniam, Ben Mueller

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
#include "MslTools.h"
#include "OptionParser.h"
#include "RegEx.h"
#include "PhiPsiStatistics.h"
#include "PhiPsiReader.h"
#include "SysEnv.h"
#include "PrositeReader.h"
#include "release.h"
#include "PolymerSequence.h"
#include "CharmmSystemBuilder.h"
#include "HydrogenBondBuilder.h"
#include "BaselineEnergyBuilder.h"
#include "PDBFragments.h"
#include "MslOut.h"
#include "CharmmEnergyCalculator.h"
#include "System.h"
#include "Residue.h"
#include "PyMolVisualization.h"
#include "AtomSelection.h"


using namespace std;
using namespace MSL;
#include <fstream>
#include <queue>

#include "aaScan.h"

struct aaDat {
  int count;
  string aa;

  aaDat(){
    count = 0;
    aa = "";
  }
  const aaDat& operator=(const aaDat &_rhs){
    if (this == &_rhs) return *this;

    count  = _rhs.count;
    aa    = _rhs.aa;
    
    return (*this);
  }

};
bool operator<( const aaDat& p1, const aaDat& p2 ) {  
      return p1.count > p2.count;
}

priority_queue<aaDat> orderMap(map<string,int> &_aMap);

static SysEnv SYSENV;
static MslOut MSLOUT("aa_scan");

int main(int argc, char *argv[]) {

	Options opt = setupOptions(argc, argv);

	// Load a PDB fragment database
	PDBFragments frags(opt.fragdb);
	frags.loadFragmentDatabase();
	//frags.setPdbDir(opt.pdbDir);
	//frags.setIncludeFullFile(opt.includeFullFile);
	//frags.printMe();
	
	vector<string> lines;
	if (MslTools::pathExtension(opt.pdb).compare("pdb") == 0){
		cout << "JUST A PDB!\n";
		lines.push_back(opt.pdb);
	} else {
		cout << "PDBLIST IS: ." << MslTools::pathExtension(opt.pdb) << "." << endl;
		MslTools::readTextFile(lines, opt.pdb);
	}

	for (uint l = 0; l < lines.size(); l++){

		// Read in query pdb
		System sys;
		sys.readPdb(lines[l]);

		// Select positions
		if (opt.sel == ""){
			opt.sel = MslTools::stringf("name CA");
		} else {
			opt.sel = MslTools::stringf("name CA and %s", opt.sel.c_str());
		}

		AtomSelection select(sys.getAtomPointers());
		AtomPointerVector ats = select.select(opt.sel);
		vector<string> positions;
		for (uint a = 0; a < ats.size(); a++){
			positions.push_back(ats[a]->getPositionId());
		}

		MSLOUT.stream() << "Number of positions to check: " << positions.size() << endl;

		// Go through each fragment length 5
		int frag_length = 5;
		for (uint p = 0; p < positions.size(); p++){

			Position &pos = sys.getPosition(positions[p]); 
			int pos_index = pos.getIndexInChain();
			Chain *ch = pos.getParentChain();

			if (pos_index - 2 + frag_length >= ch->positionSize()) continue;
			if (pos_index < 2) continue;

			MSLOUT.fprintf(stdout, "\tWORKING ON %8s-%8s-%8s\n", ch->getPosition(pos_index-2).getPositionId().c_str(), pos.getPositionId().c_str(), ch->getPosition(pos_index-2+frag_length-1).getPositionId().c_str());

			string regex = opt.regex;
			string id1 = ch->getPosition(pos_index-2).getPositionId();
			string id2 = ch->getPosition(pos_index-2+frag_length-1).getPositionId();
			int matches = frags.searchForMatchingFragmentsLinear(sys, id1, id2, regex, opt.rmsd, opt.maxMatches);

			// print out no matches for this segment?
			if (matches == 0){
				continue;
			}

			MSLOUT.fprintf(stdout, "\tNum matches: %d\n", matches);

			// matched sequences has had different formats, for now each PosId (chain-resi-icode) is a key to a sequence
			// the sequence returned is a string with all the single letter AAs that matched that position
			map<string, string> seqs = frags.getMatchedSequences();
			
			
			// Only count middle position (pos)
			//string key = MslTools::stringf("%s-%d-%s", bbAts[x]->getChainId().c_str(), bbAts[x]->getResidueNumber(), bbAts[x]->getResidueIcode().c_str());
			string key= MslTools::stringf("%s-%d-%s", pos.getChainId().c_str(), pos.getResidueNumber(), pos.getResidueIcode().c_str());
			map<string, string>::iterator it = seqs.find(key);
			
			if (it == seqs.end()){
				MSLOUT.stream() << "NO MATCHES FOR: " << key << endl;
				continue;
			}
			

			// Total number amino acids match for this position is the length of the string of single AA matches
			int total = it->second.length();

			// Code to count all amino acids by single letter code in string it->second
			map<string, int> aas;
			for (char aa : it->second) {
				aas[string(1, aa)]++;
				
			}
			
			// Print out the counts of each AA at this position
			//for (auto aa : aas) {
			//	MSLOUT.stream() << key <<"\tAA: " << aa.first << " " << aa.second << endl;
			//}
			
			
			/*
			for (it = seqs.begin(); it != seqs.end(); it++){
				if (it->first != key){
					continue;
				}

				//counts[it->second]++;
				//total++;

				// Break down each sequence into per-position per-AA counts
				//MSLOUT.stream() << "\t\t" << it->second.substr(2, 1) << endl;
				aas[it->second.substr(2, 1)]++;
			}
			*/
			MSLOUT.fprintf(stdout, "\tHERE: %d\n", total);

			// Reorder aa's at this position based on observed counts.
			priority_queue<aaDat> orderedMap = orderMap(aas);

			string outline = MslTools::stringf("%1s%d %1s %6d\t", ch->getChainId().c_str(), pos.getResidueNumber(), MslTools::getOneLetterCode(pos.getCurrentIdentity().getResidueName()).c_str(), total);
			string row = "";
			string row_counts = "";
			while (orderedMap.size() > 0){
				aaDat d = orderedMap.top();
				orderedMap.pop();

				double freq = (double)(d.count) / (double)(total) * (double)100;
				row = MslTools::stringf("  %1s %5.2f", d.aa.c_str(), freq) + row;
				row_counts = MslTools::stringf(",%1s %d %d", d.aa.c_str(), d.count, total) + row_counts;
			}
			cout << "DATA: " << outline << row << endl;
			cout << "COUNTS: " << outline << row_counts << endl;
			
			if (opt.dump != 0){

				vector<AtomContainer *> results = frags.getAtomContainers();
				vector<string> pdbnames = frags.getPDBNames();

				string outpdb = MslTools::stringf("%s_%1s%04d%s", opt.outpdb.c_str(), pos.getChainId().c_str(), pos.getResidueNumber(), pos.getResidueIcode().c_str());

				if (results.size() < opt.dump){
					opt.dump = results.size();
				}
				if (pdbnames.size() != results.size()){
					for (uint i = 0; i < results.size(); i++){
						pdbnames.push_back(MslTools::stringf("results_%06d", i));
					}
				}
				for (uint i = 0; i < opt.dump; i++){

					cout << "Writing: " << MslTools::stringf("%s_%06d.pdb", outpdb.c_str(), i) << " " << results.size() << endl;

					System newSys;
					cout << "HERE1" << endl;
					newSys.addAtoms(results[i]->getAtomPointers());
					cout << "HERE2" << endl;
					cout << "HERE2b" << pdbnames[i] << endl;
					newSys.writePdb(MslTools::stringf("%s_match_%06d_%s.pdb", outpdb.c_str(), i, pdbnames[i].c_str()));
					cout << "HERE3" << endl;
				}
			}
		}
	}
}


priority_queue<aaDat> orderMap(map<string,int> &_aMap){

  priority_queue<aaDat> results;
  map<string,int>::iterator it;
  for (it = _aMap.begin();it != _aMap.end();it++){
    aaDat a;
    a.count = it->second;
    a.aa   = it->first;
    //    MSLOUT.stream() << "Counts for "<<a.aa<<" is "<<a.count<<endl;
    results.push(a);
  }

  return results;
}


Options setupOptions(int theArgc, char * theArgv[]){
	// Create the options
	Options opt;

	// Parse the options
	OptionParser OP;
	OP.readArgv(theArgc, theArgv);
	OP.setRequired(opt.required);	
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option


	if (OP.countOptions() == 0){
		cout << "Usage: aaScan " << endl;
		cout << endl;
		cout << "\n";
		cout << "pdb PDB\n";
		cout << "fragdb DB\n";
		cout << endl;
		exit(0);
	}

	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 no pdb specified."<<endl;
		exit(1111);
	}
	opt.fragdb = OP.getString("fragdb");
	if (OP.fail()){
		cerr << "ERROR 1111 no fragdb specified."<<endl;
		exit(1111);
	}
	opt.regex = OP.getString("regex");
	if (OP.fail()){
	  opt.regex = "";
	}
	opt.sel = OP.getString("sel");
	if (OP.fail()){
	  opt.sel = "";
	}

	opt.rmsd = OP.getDouble("rmsd");
	if (OP.fail()){
	  opt.rmsd = 0.5;
	}

	// opt.dump is the number of models to print out
	opt.dump = OP.getInt("dump");
	if (OP.fail()){
	  opt.dump  = 0;
	}

	opt.outpdb = OP.getString("outpdb");
	if (OP.fail()){
	  opt.outpdb = "aaScan";
	}

	opt.pdbDir = OP.getString("pdbDir");
	if (OP.fail()){
	  opt.pdbDir = "";
	}

	opt.includeFullFile  = OP.getBool("includeAllAtoms");
	if (OP.fail()){
	  opt.includeFullFile = false;
	}
	opt.maxMatches = OP.getInt("maxMatches");
	if (OP.fail()){
	  opt.maxMatches = -1;
	}
	
	return opt;
}
