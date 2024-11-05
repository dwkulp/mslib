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

#include "RegEx.h"
#include "PolymerSequence.h"
#include "MslOut.h"

using namespace MSL;
using namespace std;

// MslOut 
static MslOut MSLOUT("RegEx");

RegEx::RegEx(){
  stype = PrimarySequence;
}

RegEx::RegEx(const RegEx &_regex){
	copy(_regex);
}

RegEx::~RegEx(){
}

void RegEx::operator=(const RegEx & _regex) {
    copy(_regex);
}

void RegEx::copy(const RegEx &_regex){
}

vector<pair<int,int> > RegEx::getResidueRanges(AtomPointerVector &_ats, string _regex){
	vector<pair<int,int> > results;
	
	string searchMe = "";
	switch (stype){
	    case PrimarySequence: 	// Atoms to 1-AA string...
	      searchMe = PolymerSequence::toOneLetterCode(_ats);
	      break;
	    case SegID:
	      for (uint i = 0; i < _ats.size();i++){
		if (_ats[i]->getSegID().length() > 0){
		  searchMe += _ats[i]->getSegID().substr(0,1);
		}
	      }
	      break;
	}

	

	//MSLOUT.stream() << "SEARCH STRING: "<<searchMe<<" ; REGEX: "<<_regex<<endl;
	cout << "SEARCH STRING: "<<searchMe<<" ; REGEX: "<<_regex<<endl;

	

	// Iterative search storing indices.
	boost::regex expression(_regex);

	// REGEX MATCH APPROACH 
	// REGEX needs to have this format ".*(WHOLE EXPRESSION).*"
	boost::smatch what;
	if (boost::regex_match(searchMe, what, expression, boost::match_extra)){

	  if (what.size() == 2){
	    pair<int,int> a;
	    a.first = what.position(1);
	    a.second = what.position(1)+what.length(1)-1;
	    results.push_back(a);
	  } else {

	    unsigned i;
	    std::cout << "** Match found **\n   Sub-Expressions:\n";
	    for(i = 2; i < what.size(); ++i) {
	      std::cout << "      $" << i << " = \"" << what[i] << " "<< what.position(i)<<"\n";
	      pair<int,int> a;
	      a.first = what.position(i);
	      a.second = what.position(i)+what.length(i)-1;
	      results.push_back(a);
	      cout << "results: "<<results.size()<<endl;
	    }
	  }
	    
	} else {
	  cout << "NO MATCH FOUND"<<endl;
	}

	// ITERATOR APPROACH
	/*
	//boost::sregex_token_iterator r1(searchMe.begin(),searchMe.end(),expression,-1);
	boost::sregex_iterator r1(searchMe.begin(),searchMe.end(),expression);
	boost::sregex_iterator r2;
	if (r1 == r2){
	  cout << "NO MATCH\n";
	}
	while (r1 != r2){
		cout << "MATCH: "<<*r1<<" "<<(*r1).position()<<" "<<(*r1).length()<<endl;

		pair<int,int> a;
		a.first = (*r1).position();
		a.second = (*r1).position()+(*r1).length()-1;
		results.push_back(a);

		*r1++;
	}
	*/
	return results;

}

vector<pair<int,int> > RegEx::getResidueRanges(Chain &_ch, string _regex){

	vector<pair<int,int> > results;
	
	
	string searchMe = "";
	switch (stype){
	    case PrimarySequence: 	// Chain to 1-AA string...
	      searchMe = PolymerSequence::toOneLetterCode(_ch.getAtomPointers());
	      break;
	    case SegID:
	      for (uint i = 0; i < _ch.positionSize();i++){
		Position &pos = _ch.getPosition(i);
		if (pos.atomExists("CA") && pos.getAtom("CA").getSegID().length() > 0){
		  searchMe += pos.getAtom("CA").getSegID().substr(0,1);
		}
	      }
	      break;
	}

	

	//MSLOUT.stream() << "SEARCH STRING: "<<searchMe<<" ; REGEX: "<<_regex<<endl;
	cout << "SEARCH STRING: "<<searchMe<<" ; REGEX: "<<_regex<<endl;
	// Iterative search storing indices.
	boost::regex expression(_regex);

	//boost::sregex_token_iterator r1(searchMe.begin(),searchMe.end(),expression,-1);
	boost::sregex_iterator r1(searchMe.begin(),searchMe.end(),expression);
	boost::sregex_iterator r2;
	if (r1 == r2){
	  cout << "NO MATCH\n";
	}
	/*
	MSLOUT.stream() << "r1: "<<(*r1).size()<<endl;
	for (uint i =0 ;i < (*r1).size();i++){
		MSLOUT.stream() << "M: "<<(*r1).position()<<endl;
		*r1++;
	}
	*/
	while (r1 != r2){
		MSLOUT.stream() << "MATCH: "<<*r1<<" "<<(*r1).position()<<" "<<(*r1).length()<<endl;

		pair<int,int> a;
		a.first = (*r1).position();
		a.second = (*r1).position()+(*r1).length()-1;
		results.push_back(a);

		*r1++;
	}

	return results;

}
