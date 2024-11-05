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


#include "AtomSelection.h"

using namespace MSL;
using namespace std;


AtomSelection::AtomSelection(){
	data = NULL;
	debug = false;
	storedSelections.clear();
	storedSelections["_EMPTY_"]; // create an empty selection
	storedSelections["ALL"];
}

AtomSelection::AtomSelection(AtomPointerVector &_data){
	data = &_data;
	debug = false;
	storedSelections.clear();
	storedSelections["_EMPTY_"]; // create an empty selection
	storedSelections["ALL"] = *data; // create the ALL selection
}

AtomSelection::~AtomSelection(){
}


bool AtomSelection::selectionExists(string _selectName){
	Hash<string,AtomPointerVector>::Table::iterator it = storedSelections.find(_selectName);
	return (it != storedSelections.end());
}
AtomPointerVector& AtomSelection::getSelection(string _selectName){
	_selectName = MslTools::toUpper(_selectName);
	Hash<string,AtomPointerVector>::Table::iterator it = storedSelections.find(_selectName);
	if (it != storedSelections.end()){
		//a = &(it->second);
		return it->second;
	}

	return storedSelections["_EMPTY_"];
}

unsigned int AtomSelection::selectionSize(string _selectName){
	_selectName = MslTools::toUpper(_selectName);
	Hash<string,AtomPointerVector>::Table::iterator it = storedSelections.find(_selectName);
	if (it != storedSelections.end()){
		//a = &(it->second);
		return it->second.size();
	}

	return 0;
}



  

AtomPointerVector& AtomSelection::inverseSelect(string _selectString, bool _selectAllAtoms){

  string newSelectStr = "not ("+_selectString+")";
  return select(newSelectStr,_selectAllAtoms);
  
}

// THE FOLLOWING PRECOMPILER DIRECTIVE DECIDES WHAT TO COMPILE BETWEEN TWO IMPLEMENTATIONS OF THE SELECTION METHOD
#ifndef __TESTING__	
// DEFAULT
// ALESSANDRO'S CODE
AtomPointerVector& AtomSelection::select(string _selectString, bool _selectAllAtoms){
	
	// A selection name may be given, separated from the selection string by a comma ("res37CA, resi 37 AND name CA")
	_selectString = MslTools::trim(_selectString);
	vector<string> selectToks = MslTools::tokenizeAndTrim(_selectString, ",");

	string name = "";
	string selectString = selectToks[0];
	if (selectToks.size() > 1){
		// a selection name was given
		name = MslTools::toUpper(selectToks[0]);
		selectString = selectToks[1];
	}
	if (debug) {
		cout << "Selection string, name: \"" << name << "\", logic: \"" << selectString << "\"" << endl;
	}

	// this statement clears the selection or it creates the new selection if it didn't exist (and emtpy AtomPointerVector)
	storedSelections[name].clear();

	if (data == NULL || data->size() == 0) {
		// Nothing to do...
		return storedSelections[name];
	}

	// Chcek if a WITHIN X OF was given
	size_t wpos = selectString.find("WITHIN",0);
	size_t opos = selectString.find("OF",0);
	if (wpos != std::string::npos && opos != std::string::npos) {
		/*********************************************************
		 *  It is a WITHIN X OF selection (i.e. NAME CA WITHIN 5 OF NAME CB)
		 *
		 *  1) split the statement, getting the selections strings before (NAME CA)
		 *     and after (NAME CB)
		 *  2) create a temporary selection for the second string (NAME CB)
		 *  3) find all atoms that are within X Angstrom of the temporary
		 *     selection (including those in the temporary selection)
		 *
		 *********************************************************/
		if (wpos > opos) {
			cerr << "WARNING 27473: selection statement with WITHIN <A> OF " << _selectString << " not valid in AtomPointerVector& AtomSelection::select(string _selectString, bool _selectAllAtoms)" << endl;
			return storedSelections[name];
		}
		string seleString1 = selectString.substr(0, wpos-1);
		string withinPart = selectString.substr(wpos, opos-wpos+2);
		string seleString2 = selectString.substr(opos+3);

		vector<string> withinTokens = MslTools::tokenizeAndTrim(withinPart);
		if (withinTokens.size() != 3 || withinTokens[0] != "WITHIN" || withinTokens[2] != "OF") {
			cerr << "WARNING 27478: selection statement with WITHIN <A> OF " << _selectString << " not valid in AtomPointerVector& AtomSelection::select(string _selectString, bool _selectAllAtoms)" << endl;
			return storedSelections[name];
		}
		double radius = MslTools::toDouble(withinTokens[1]);

		// create the temporary selection from string #2
		AtomPointerVector sele2atoms = logicalSelect(seleString2, "_TMP1_", *data, _selectAllAtoms);

		// find any atom that is within the range of selection 2
		AtomPointerVector aroundSele2;
		for (AtomPointerVector::iterator avIt = data->begin();avIt != data->end();avIt++){
			if (!(*avIt)->hasCoor()) {
				// skip atoms that do not have coordinates
				continue;
			}
			if((*avIt)->getSelectionFlag("_TMP1_")) {
				// it is part of the sele1, it goes in
				aroundSele2.push_back(*avIt);
				continue;
			}
			for (AtomPointerVector::iterator avJt = sele2atoms.begin();avJt != sele2atoms.end();avJt++){
				if (!(*avJt)->hasCoor()) { 
					// skip atoms that do not have coordinates
					continue;
				}
				if ((*avIt)->distance(**avJt) < radius) {
					// atom within the radius of sele 2
					aroundSele2.push_back(*avIt);
					break;
				}
			}
		}

		// delete the temporary selection and return
		clearStoredSelection("_TMP1_");
		return logicalSelect(seleString1, name, aroundSele2, _selectAllAtoms);

	} else {
		// regular logical statement, parse, select and return
		return logicalSelect(selectString, name, *data, _selectAllAtoms);
	}
}

AtomPointerVector& AtomSelection::logicalSelect(string _selectString, string _name, AtomPointerVector & _atoms, bool _selectAllAtoms){

	if (!cond.setLogic(_selectString)) {
		cerr << "WARNING 32773: selection statement " << _selectString << " not valid in AtomPointerVector& AtomSelection::logicalSelect(string _selectString, string _name, AtomPointerVector & _atoms, bool _selectAllAtoms)" << endl;
		return storedSelections[_name];
	}
	if (debug){
		cout << "Reconstructed selection logic: " << cond.printLogicalConditions() << endl;
	}

	for (AtomPointerVector::iterator avIt = _atoms.begin();avIt != _atoms.end();avIt++){
		//  let's reset the query status of the condition
		cond.restartQuery();

		// first check if the atoms should be discarded because inactive
		if (_selectAllAtoms || (*avIt)->getActive()) {

			bool success = true;
			// answer all the question asked by the LogicalCondition
			while (!cond.logicComplete()) {
				vector<string> tokens = cond.getLogicalCondition();
				if (tokens.size() < 1) {
					// invalid condition
					cerr << "WARNING 32778: blank condition (ignored) in selection statement " << _selectString << " in AtomPointerVector& AtomSelection::logicalSelect(string _selectString, string _name, AtomPointerVector & _atoms, bool _selectAllAtoms)" << endl;
					continue;
				}

				// check if the atom satisfies the condition
				bool condition = false;
				if (tokens[0] == "ALL" && tokens.size() == 1) {
					// all atoms
					condition = true;
				} else if (tokens[0] == "NAME" && tokens.size() == 2) {
					vector<string> vals=MslTools::tokenize(tokens[1], "+"); // split if it is multi CA+CB+CG
					for (unsigned int i=0; i<vals.size(); i++) {
						if (vals[i] == (*avIt)->getName()) {
							condition = true;
							break;
						} 
					}
				} else if (tokens[0] == "RESI" && tokens.size() == 2) {
					vector<string> vals=MslTools::tokenize(tokens[1], "+"); // split if it is multi 7+18+22
					for (unsigned int i=0; i<vals.size(); i++) {
						vector<string> range=MslTools::tokenize(vals[i], "-"); // split if it is ranges 
						if (range.size() == 2) {
							// range: resi 6-9
							int start = MslTools::toInt(range[0]);
							int end = MslTools::toInt(range[1]);
							if (start <= (*avIt)->getResidueNumber() && end >= (*avIt)->getResidueNumber()) {
								condition = true;
								break;
							} 
						} else {
							// note this support also a possible insertion code by using the positionId with
							// a skip-level of 1 (skip chain), i.e. "37" or "37A" (the above range doesn't)
							if (vals[i] == (*avIt)->getPositionId(1)) {
								condition = true;
								break;
							} 
						}
//<<<<<<< .mine
					}
				} else if (tokens[0] == "RESN" && tokens.size() == 2) {
					vector<string> vals=MslTools::tokenize(tokens[1], "+"); // split if it is multi ALA+LEU+VAL
					for (unsigned int i=0; i<vals.size(); i++) {
						if (vals[i] == (*avIt)->getResidueName()) {
//=======
//					} else if (tokens[0] == "CHAIN" && tokens.size() == 2) {
//						vector<string> vals=MslTools::tokenize(tokens[1], "+"); // split if it is multi A+B+C
//						for (unsigned int i=0; i<vals.size(); i++) {
//							if (vals[i] == (*avIt)->getChainId()) {
//								condition = true;
//								break;
//							} 
//						}
//					} else if (tokens[0] == "HASCRD" || tokens[0] == "HASCOOR") {
//						bool cmp = true;
//						if (tokens.size() == 2) {
//							// HASCRD FALSE or HASCRD 0 or HASCRD TRUE or HASCRD 1
//							cmp = MslTools::toBool(tokens[1]);
//						}
//						if ((*avIt)->hasCoor() == cmp) {
//>>>>>>> .r450
							condition = true;
							break;
						} 
					}
				} else if (tokens[0] == "CHAIN" && tokens.size() == 2) {
					vector<string> vals=MslTools::tokenize(tokens[1], "+"); // split if it is multi A+B+C
					for (unsigned int i=0; i<vals.size(); i++) {
						if (vals[i] == (*avIt)->getChainId()) {
							condition = true;
							break;
						} 
					}
				} else if (tokens[0] == "HASCRD" || tokens[0] == "HASCOOR") {
					bool cmp = true;
					if (tokens.size() == 2) {
						// HASCRD FALSE or HASCRD 0 or HASCRD TRUE or HASCRD 1
						cmp = MslTools::toBool(tokens[1]);
					}
					if ((*avIt)->hasCoor() == cmp) {
						condition = true;
					}
				} else if (tokens[0] == "CHARGE") {
					int cmp = 0;
					if (tokens.size() == 2) {
						cmp = MslTools::toInt(tokens[1]);
						if ((*avIt)->getCharge() == cmp) {
						  condition = true;
						}
					}
				} else if (tokens[0] == "TEMPFACTOR") {
					int cmp = 0;
					if (tokens.size() == 2) {
						cmp = MslTools::toInt(tokens[1]);
						if ((*avIt)->getTempFactor() == cmp) {
						  condition = true;
						}
					}
				} else if (tokens.size() == 1 && storedSelections.find(tokens[0]) != storedSelections.end()) {
					// another selection
					if ((*avIt)->getSelectionFlag(tokens[0])) {
						condition = true;
					}
				} else {
					cerr << "WARNING 32783: unrecognized condition (ignored) in selection statement " << _selectString << " in AtomPointerVector& AtomSelection::logicalSelect(string _selectString, string _name, AtomPointerVector & _atoms, bool _selectAllAtoms)" << endl;
					success = false;
					break;
				}
//<<<<<<< .mine
				//if (debug) {
				//	cout << "CONDITION: " << MslTools::joinLines(tokens) << " " << condition << endl;
				//}
				// pass the info to the LogicalCondition, the object keeps track of the logic
				cond.setLogicalConditionValue(condition);
//=======
//				// time to decide if this atoms is selected or not
//				bool selected = cond.getOverallBooleanState();
//				if (debug) {
//					cout << cond.printLogicalConditionsWithValues() << endl;
//				}
//				if (selected) {
//					storedSelections[name].push_back(*avIt);
//				}
//				// turn the atom's selection flag on/off
//				(*avIt)->setSelectionFlag(name, selected);
//
//>>>>>>> .r450
			}
			if (!success) {
				// the atom selection syntax was problematic
				break;
			}
			// time to decide if this atoms is selected or not
			bool selected = cond.getOverallBooleanState();
			if (debug) {
				cout << cond.printLogicalConditionsWithValues() << endl;
			}
			if (selected) {
				storedSelections[_name].push_back(*avIt);
			}
			// turn the atom's selection flag on/off
			(*avIt)->setSelectionFlag(_name, selected);

		}

	}


	return storedSelections[_name];

}
// END OF ALESSANDRO'S CODE

#else
// IF $MSL_TESTING IS TRUE
// DAN'S CODE
AtomPointerVector& AtomSelection::select(string _selectString, bool _selectAllAtoms){
	
	// Tokenize select string based on ',' character 
	vector<string> selectToks = MslTools::tokenize(_selectString, ",");

	string name         = MslTools::toUpper(selectToks[0]);
	string selectString = selectToks[0];
	if (selectToks.size() > 1){
		selectString = selectToks[1];
	}

	//AtomPointerVector *a = NULL;

	// Check for no data
	if (data == NULL) {
		return storedSelections["_EMPTY_"];
		//cerr << "WARNING no AtomPointerVector to select from, most likely you forgot to give this AtomSelection object a starting AtomPointerVector in its constructor."<<endl;
		//exit(-1);
	}


	// Create AtomPointerVector for data
	AtomPointerVector tmp;

	/*
	  Check for complex selections:  
	  1.   SEL1 within 5 of SEL2 
	*/


	int wpos = selectString.find("WITHIN",0);
	int opos = selectString.find("OF",0);
	if (wpos != std::string::npos && opos != std::string::npos) {

		// Complex Selection form 'SEL1 WITHIN X of SEL2'
		//   Create 2 logical parsers.  Create two tmp AtomPointerVectors.  Do combination.
			
		string sel1 = MslTools::trim(selectString.substr(0,wpos));
		double dist = MslTools::toDouble(MslTools::trim(selectString.substr(wpos+6,wpos-opos)),"AtomSelection::select is trying to find a distance in selection statement");
		string sel2 = MslTools::trim(selectString.substr(opos+2));

		//cout << "sel1: '"<<sel1<<"'"<<endl;
		//cout << "dist: '"<<dist<<"'"<<endl;
		//cout << "sel2: '"<<sel2<<"'"<<endl;

		// Parse the Logic
		LogicalParser lp1;
		lp1.setDebugFlag(debug);
		lp1.setLogicStatementInFix(sel1);
		lp1.parse();

		LogicalParser lp2;
		lp2.setDebugFlag(debug);
		lp2.setLogicStatementInFix(sel2);
		lp2.parse();


		// Locally store sel1,sel2
		AtomPointerVector tmp1;
		AtomPointerVector tmp2;
		AtomPointerVector::iterator avIt;
		for (avIt = data->begin();avIt != data->end();avIt++){

			if (lp1.eval(**(avIt)) && (_selectAllAtoms || (**(avIt)).getActive())){
				tmp1.push_back(*avIt);
			}

			if (lp2.eval(**(avIt)) && (_selectAllAtoms || (**(avIt)).getActive())) {
				tmp2.push_back(*avIt);
			}

		}

		AtomPointerVector::iterator avIti;
		AtomPointerVector::iterator avItj;

		for (avIti = tmp1.begin(); avIti != tmp1.end(); avIti++){

			bool withinDistance = false;
			for (avItj = tmp2.begin(); avItj != tmp2.end(); avItj++){

				if ((*avIti)->distance((**avItj)) <= dist) {
					withinDistance = true;
					break;
				}
			}

			if (withinDistance){
				(*avIti)->setSelectionFlag(name, true);
				tmp.push_back(*avIti);
			} else {
				(*avIti)->setSelectionFlag(name, false);
			}
				
		}

	} else {


		// Parse the Logic
		if (debug){
			cout << "Selection Statement: "<<selectString<<endl;
		}

		LogicalParser lp;
		lp.setDebugFlag(debug);
		lp.setLogicStatementInFix(selectString);
		lp.parse();

		if (debug){
			lp.printLogicTree();
		}


		AtomPointerVector::iterator avIt;
		for (avIt = data->begin();avIt != data->end();avIt++){
			//cout  << "Data: "<<**avIt<<endl;

			if (lp.eval(**(avIt)) && (_selectAllAtoms || (**(avIt)).getActive())) {
				(*avIt)->setSelectionFlag(name, true);
				tmp.push_back(*avIt);
			} else {
				(*avIt)->setSelectionFlag(name, false);
			}

		}
	}


	storedSelections[name] = tmp; // this does a copy
	//a =  &storedSelections[name];
		



	// BUG!!! Returning a local pointer by reference!!!
	//return (*a);
	return storedSelections[name];
}

#endif


