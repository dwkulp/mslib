/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2012 The MSL Developer Group (see README.TXT)
 MSL Libraries: http://msl-libraries.org

If used in a scientific publication, please cite: 
 Kulp DW, Subramaniam S, Donald JE, Hannigan BT, Mueller BK, Grigoryan G and 
 Senes A "Structural informatics, modeling and design with a open source 
 Molecular Software Library (MSL)" (2012) J. Comput. Chem, 33, 1645-61 
 DOI: 10.1002/jcc.22968

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

#ifndef PDBFRAGMENTS_H
#define PDBFRAGMENTS_H

#include <string>

#include "AtomPointerVector.h"
#include "AtomContainer.h"
#include "System.h"

namespace MSL { 
class PDBFragments{

	public:
		PDBFragments();
		PDBFragments(std::string _fragDbFile, std::string _BBQTableForBackboneAtoms="");
		~PDBFragments();


		void loadFragmentDatabase();
		void setFragDB(std::string _fragdb);
		void setPdbDir(std::string _pdbdir);
		void setBBQTable(std::string _table);

        int searchForMatchingFragmentsLinear(System &_sys, string &_startRes, string &_endRes, string _regex, double _rmsdTol,int _maxFrags=-1);
		int searchForMatchingFragmentsStems(System &_sys, std::vector<std::string> &_stemResidues,int _numResiduesInFragment=-1,std::string _regex="",double _rmsdTol=0.5);
		int searchForMatchingFragmentsSpots(System &_sys, std::vector<std::string> &_stemResidues,int _maxResiduesInFragment, double _rmsdTol);

		// Function searches for fragments within "_searchSys" rather than the included fragDB.
		int searchForMatchingDualFragments(System &_sys1, std::vector<std::string> &_stemResidues1,
		System &_sys2, std::vector<std::string> &_stemResidues2,
		System &_searchSys, int _loop1min, int _loop1max, int _loop2min, int _loop2max, double _distanceStem1=0.0, double _distanceStem2=0.0,
		double _stemRmsdTol=0.3,double _totalRmsdTol=0.6,bool _matchFirstStemOnly=false);

		// Function searches for fragments inside fragDB.
		int searchForMatchingDualFragments(System &_sys1, std::vector<std::string> &_stemResidues1,
		System &_sys2, std::vector<std::string> &_stemResidues2,
		int _loop1min, int _loop1max, int _loop2min, int _loop2max, double _distanceStem1=0.0, double _distanceStem2=0.0,
		double _stemRmsdTol=0.3,double _totalRmsdTol=0.6,bool _matchFirstStemOnly=false);

		vector<AtomContainer *> & getAtomContainers();
		AtomPointerVector getAtomPointers();
		map<std::string,std::string> & getMatchedSequences();

		enum dbAtoms { caOnly=0, allAtoms=1 };
		enum searchType { linear=0, stemOnly=1, discreteSpots=2 };
		void printMe();

		void setIncludeFullFile(bool _flag);

		void setMaxFragments(int _maxNumFrags);
		int getMaxFragments();

		void setStoreResfileLines(bool _flag) { storeResfileLinesFlag = _flag; }
		vector<string> getResfileLines() { return resfile_lines; }
		vector<string> getPDBNames() { return pdbnames; }
	private:
		std::string fragDbFile;
		dbAtoms fragType;
		std::string bbqTable;
		AtomPointerVector fragDB;
		std::string pdbDir;
		bool includeFullFile;
		int maxNumFrags;
		bool storeResfileLinesFlag;
		vector<string> resfile_lines;
		vector<string> pdbnames;
		map<std::string,std::string> matchedSequences;
		vector<AtomContainer *> lastResults;
};

inline void PDBFragments::setFragDB(std::string _fragdb) { fragDbFile = _fragdb;}
inline void PDBFragments::setPdbDir(std::string _pdbdir) { pdbDir = _pdbdir;}
inline void PDBFragments::setBBQTable(std::string _table) { bbqTable = _table;}
inline vector<AtomContainer*> & PDBFragments::getAtomContainers() { return lastResults;}		
inline AtomPointerVector PDBFragments::getAtomPointers() { 
  AtomPointerVector ats;
  for (uint i =0;  i < lastResults.size();i++){
    ats += lastResults[i]->getAtomPointers();
  }
  return ats;
}
inline PDBFragments::PDBFragments() { 	fragType   = caOnly; pdbDir = ""; fragDbFile = ""; bbqTable="";}
inline PDBFragments::PDBFragments(std::string _fragDbFile,std::string _BBQTableForBackboneAtoms) {
	fragDbFile = _fragDbFile;
	pdbDir = "";

	if (_BBQTableForBackboneAtoms == ""){
		fragType   = allAtoms;
	} else {
		fragType   = caOnly;
	}
	bbqTable = _BBQTableForBackboneAtoms;

}
inline PDBFragments::~PDBFragments() {
}
inline void PDBFragments::loadFragmentDatabase(){
	fragDB.load_checkpoint(fragDbFile);
	cout << "FragDB: "<<fragDB.getName()<<" has "<<fragDB.size()<<" atoms."<<endl;
	if (fragDB.getName() == "ca-only"){
	  fragType = caOnly;
	}

	if (fragDB.getName() == "allatom"){
	  fragType = allAtoms;
	}

	// Reset the active conformation as this will reset the currentIterator for atom pointer vectors
	fragDB.setActiveConf(0);

}

inline map<std::string,std::string> & PDBFragments::getMatchedSequences(){
  return matchedSequences;
}
inline void PDBFragments::printMe(){
	for (uint i = 0; i < fragDB.size();i++){
		std::cout << fragDB(i).toString()<<fragDB(i).getSegID()<<std::endl;
	}
}

inline void PDBFragments::setIncludeFullFile(bool _flag){
  includeFullFile = _flag;
}
 
inline void PDBFragments::setMaxFragments(int _maxNumFrags){ maxNumFrags = _maxNumFrags;}
inline int PDBFragments::getMaxFragments(){ return maxNumFrags;}
}

#endif
