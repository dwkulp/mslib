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

#include "Residue.h"
#include "Position.h"
#include "System.h"

Residue::Residue() : Selectable<Residue>(this) {
	setup("DUM", 1, "", "A");
}

Residue::Residue(string _resName, int _resNum, string _icode) : Selectable<Residue>(this){
	setup(_resName, _resNum, _icode, "A");
}

Residue::Residue(const AtomVector & _atoms, string _resName, int _resNum, string _icode) : Selectable<Residue>(this){
	setup(_resName, _resNum, _icode, "A");
	addAtoms(_atoms);
}

Residue::Residue(const Residue & _residue) : Selectable<Residue>(this){
	pParentPosition = NULL;
	copy(_residue);

	addSelectableFunctions();

}

Residue::~Residue() {
	deletePointers();
}

void Residue::addSelectableFunctions(){

	addStringFunction("RESN", &Residue::getResidueName);
	addIntFunction("RESI", &Residue::getResidueNumber);
	addStringFunction("ICODE", &Residue::getResidueIcode);
	addStringFunction("CHAIN", &Residue::getChainId);


	//addBoolFunction("CA", &Residue::exists("CA"));

}

void Residue::operator=(const Residue & _residue) {
	copy(_residue);
}

void Residue::setup(string _resName, int _resNum, string _insertionCode, string _chainId) {
	pParentPosition = NULL;
	residueName = _resName;
	residueNumber = _resNum;
	residueIcode = _insertionCode;
	chainId = _chainId;

	addSelectableFunctions();

}

void Residue::copy(const Residue & _residue) {
	deletePointers();
	//residueName = _residue.residueName;
	setResidueName(_residue.residueName); // this updates the Identity map
	residueNumber = _residue.residueNumber;
	residueIcode = _residue.residueIcode;
	chainId = _residue.getChainId();
	for (vector<AtomGroup*>::const_iterator k=_residue.electrostaticGroups.begin(); k!=_residue.electrostaticGroups.end(); k++) {
		// for each group
		for (AtomGroup::const_iterator l=(*k)->begin(); l!=(*k)->end(); l++) {
			// k-electrostaticGroups.begin() (iterator subtraction) is the integer index of the group
			addAtom(**l);
		}
	}
}

void Residue::deletePointers() {
	for (AtomVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		delete *k;
	}
	atoms.clear();
	for (vector<AtomGroup*>::iterator k=electrostaticGroups.begin(); k!=electrostaticGroups.end(); k++) {
		delete *k;
	}
	electrostaticGroups.clear();
	atomMap.clear();
}

void Residue::addAtom(const Atom & _atom) {
	/*
	CartesianPoint c(_atom.getX(), _atom.getY(), _atom.getZ());
	unsigned int group = _atom.getGroupNumber();
	addAtom(_atom.getName(), c, group);
	*/
	string name = _atom.getName();
	map<string, Atom*>::iterator foundAtom=atomMap.find(name);
	unsigned int group = _atom.getGroupNumber();

	if (foundAtom==atomMap.end()) {
		/******************************************
		 *  An atom with _atom's name DOES NOT exist: 
		 *   - add another atom
		 ******************************************/
		atoms.push_back(new Atom(_atom));
		while (group >= electrostaticGroups.size()) {
			// add groups if _group is more then the last index
			electrostaticGroups.push_back(new AtomGroup(this));
		}
		// add the atom to the atom to the electrostatic group
		// and to the map
		(electrostaticGroups[group])->push_back(atoms.back());
		atomMap[name] = atoms.back();



	} else {
		/******************************************
		 *  An atom with _atom's name DOES ALREADY exist: 
		 *   - add alternate coordinates to the atom
		 *  
		 *  Note: in this case we are disregarding the _group
		 ******************************************/
		CartesianPoint c(_atom.getX(), _atom.getY(), _atom.getZ());
		(*foundAtom).second->addAltConformation(c);
	}


}

void Residue::addAtom(string _name, const CartesianPoint & _coor, size_t _group) {
	//cout << "Adding atom " << _name << " w/ cartesian point to group " << _group << endl;

	map<string, Atom*>::iterator foundAtom=atomMap.find(_name);

	if (foundAtom==atomMap.end()) {
		/******************************************
		 *  An atom with _atom's name DOES NOT exist: 
		 *   - add another atom
		 ******************************************/
		atoms.push_back(new Atom(_name, _coor));
		while (_group >= electrostaticGroups.size()) {
			// add groups if _group is more then the last index
			electrostaticGroups.push_back(new AtomGroup(this));
		}
		// add the atom to the atom to the electrostatic group
		// and to the map
		//cout << "UUU Atom add atom " << atoms.back()->getName() << " atoms" << endl;
		//cout << "UUU _group = " << _group << " electrostaticGroups.size() = " << electrostaticGroups.size() << endl;
		(electrostaticGroups[_group])->push_back(atoms.back());
		atomMap[_name] = atoms.back();



	} else {
		/******************************************
		 *  An atom with _atom's name DOES ALREADY exist: 
		 *   - add alternate coordinates to the atom
		 *  
		 *  Note: in this case we are disregarding the _group
		 ******************************************/
		//cout << "UUU Atom add alt conf to atom " << (*foundAtom).second->getName() << " atoms" << endl;
		(*foundAtom).second->addAltConformation(_coor);
	}
	//cout << "UUU Residue return" << endl;
}		

void Residue::addAtoms(const AtomVector & _atoms) {
	for (AtomVector::const_iterator k=_atoms.begin(); k!=_atoms.end(); k++) {
		addAtom(*(*k)); // *k is an atom pointer, *(*k) and atom
	}
}

void Residue::addAltConformationToAtom(string _name, const CartesianPoint & _coor) {
	//  we add this function for clarity but addAtoms already
	//  does the job
	addAtom(_name, _coor, 0);
}

bool Residue::removeAtom(string _name) {
	map<string, Atom*>::iterator foundAtom=atomMap.find(_name);

	if (foundAtom!=atomMap.end()) {
		// erase from its electrostatic group
		for (vector<AtomGroup*>::iterator k=electrostaticGroups.begin(); k!=electrostaticGroups.end(); k++) {
			// for each group
			for (AtomGroup::iterator l=(*k)->begin(); l!=(*k)->end(); l++) {
				// for each Atom pointer
				if ((*foundAtom).second == *l) {
					// found it, remove it from the group
					(*k)->erase(l);
				}
			}
		}
		// erase from the atoms list
		for (AtomGroup::iterator k=atoms.begin(); k!=atoms.end(); k++) {
			if ((*foundAtom).second == *k) {
				// deallocate from memory and remove from list
				delete *k;
				atoms.erase(k);
				// erase from the map
				atomMap.erase(foundAtom);
				return true;
			}
		}
	}
	return false;
}

void Residue::removeAllAtoms() {
	deletePointers();
}

void Residue::setResidueNumber(int _resNum) {
	if (pParentPosition != NULL) {
		pParentPosition->setResidueNumber(_resNum);
	} else {
		residueNumber = _resNum;
	}
}

int Residue::getResidueNumber() const {
	if (pParentPosition != NULL) {
		return pParentPosition->getResidueNumber();
	} else {
		return residueNumber;
	}
}

void Residue::setResidueIcode(string _icode) {
	if (pParentPosition != NULL) {
		pParentPosition->setResidueIcode(_icode);
	} else {
		residueIcode = _icode;
	}
}

string Residue::getResidueIcode() const {
	if (pParentPosition != NULL) {
		return pParentPosition->getResidueIcode();
	} else {
		return residueIcode;
	}
}

void Residue::setChainId(string _chainId) {
	if (pParentPosition != NULL) {
		pParentPosition->setChainId(_chainId);
	} else {
		chainId = _chainId;
	}
}

string Residue::getChainId() const {
	if (pParentPosition != NULL) {
		return pParentPosition->getChainId();
	} else {
		return chainId;
	}
}

// the name space defines if the atoms and residues are named with
// PDB or charmm names (or whatever else)
void Residue::setNameSpace(string _nameSpace) {
	if (pParentPosition != NULL) {
		pParentPosition->setNameSpace(_nameSpace);
	} else {
		nameSpace = _nameSpace;
	}
}

string Residue::getNameSpace() const {
	if (pParentPosition != NULL) {
		return pParentPosition->getNameSpace();
	} else {
		return nameSpace;
	}
}

bool Residue::getActive() const {
	if (pParentPosition == NULL) {
		return true;
	} else {
		// the residue is active is its address is the same of the
		// currently active residue in the parent Position
		return this == &(pParentPosition->getCurrentIdentity());
	}
}



void Residue::updateAtomMap(Atom * _atom) {
	// if an atom changes its name it calls this function
	// to update the map (through its group)
	for (map<string, Atom*>::iterator k=atomMap.begin(); k!=atomMap.end(); k++) {
		if (k->second == _atom) {
			if (k->first != _atom->getName()) {
				atomMap.erase(k);
				atomMap[_atom->getName()] = _atom;
			}
		}
	}
}

void Residue::updatePositionMap() {
	if (pParentPosition != NULL) {
		pParentPosition->updateResidueMap(this);
	}
}


/***************************************************
 *  Manage the alternative conformations of the atoms
 ***************************************************/
void Residue::setActiveConformation(size_t _i) {
	// set the atoms to their i-th conformation (or
	// 0-th if they do not have it
	for (AtomVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		if ((*k)->getNumberOfAltConformations() > _i) {
			(*k)->setActiveConformation(_i);
		} else {
			(*k)->setActiveConformation(0);
		}
	}
}

/*
unsigned int Residue::getActiveConformation() const {
	/ ***************************************************
	 *  A residue is in a alt conformation different from
	 *  the 0-th if all its atoms that have more than 1
	 *  alt conf are consistent in the same conformation
	 *  (those that have just one conf do not count)
	 *
	 *  NOT USED BECAUSE IT IS AMBIGUOUS: 
	 *   - if out > 0, the residue is consistent
	 *     in a alt conf
	 *   - but if out == 0, all atoms might be in the 0-th
	 *     conf or the atoms might be inconsistent
	 *************************************************** /
	unsigned int out = 0;
	bool foundAtomWithMultipleConf = false;
	for (AtomVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		if ((*k)->getNumberOfAltConformations() == 1) {
			// atoms with just 1 conformation do not count
			continue;
		}
		if (!foundAtomWithMultipleConf) {
			out = (*k)->getActiveConformation();
			foundAtomWithMultipleConf = true;
		} else {
			if ((*k)->getActiveConformation() != out) {
				// not consistent, return 0
				return 0;
			}
		}
	}
	return out;
}
*/

unsigned int Residue::getNumberOfAltConformations() const {
	// the highest number of alt conf among the atoms
	unsigned int max = 0;
	for (AtomVector::const_iterator k=atoms.begin(); k!=atoms.end(); k++) {
		if (k==atoms.begin() || (*k)->getNumberOfAltConformations() > max) {
			max = (*k)->getNumberOfAltConformations();
		}
	}
	return max;
}

void Residue::addAltConformation() {
	for (AtomVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		(*k)->addAltConformation();
	}
}

void Residue::addAltConformation(const vector<CartesianPoint> & _points) {
	for (AtomVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		if (_points.size() > k - atoms.begin()) {
			(*k)->addAltConformation(_points[k - atoms.begin()]);
		} else {
			// if we did not have enough points, default it
			(*k)->addAltConformation();
		}
	}
}

void Residue::removeAltConformation(size_t _i) {
	for (AtomVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		(*k)->removeAltConformation(_i);
	}
}

void Residue::removeAllAltConformations() {
	for (AtomVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		(*k)->removeAllAltConformations();;
	}
}


Chain * Residue::getParentChain() const {
	if (pParentPosition != NULL) {
		return pParentPosition->getParentChain();
	} else {
		return NULL;
	}
}

System * Residue::getParentSystem() const {
	if (pParentPosition != NULL) {
		return pParentPosition->getParentSystem();
	} else {
		return NULL;
	}
}

unsigned int Residue::getIdentityIndex() {
	if (pParentPosition != NULL) {
		return pParentPosition->getIdentityIndex(this);
	} else {
		return 0;
	}
}

vector<int> Residue::findNeighbors(double _distance){

	System *p = getParentSystem();

	if (p == NULL){
		cerr << "ERROR 1966 Residue::findNeighbors() has a NULL parent System.\n";
		exit(1966);
	}


	vector<int> result;
	for (uint i = 0 ; i < p->residueSize();i++){

		Residue &r = p->getResidue(i);

		// Skip over this residue
		if (&r == this){
			continue;
		}

		if (distance(r,"CENTROID") < _distance){
			result.push_back(i);
		}
		
	}

	return result;
}
