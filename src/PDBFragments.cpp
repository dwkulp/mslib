/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2013 The MSL Developer Group (see README.TXT)
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


#include "PDBFragments.h"
#include "BBQTable.h"
#include "Transforms.h"
#include "AtomSelection.h"
#include "MslExceptions.h"
#include "MslOut.h"

// BOOST Includes
#include <boost/regex.hpp>


static MslOut MSLOUT("PDBFragments");
using namespace MSL;
using namespace std;

int PDBFragments::searchForMatchingDualFragments(System &_sys1, std::vector<std::string> &_stemResidues1,
						 System &_sys2, std::vector<std::string> &_stemResidues2,
						 int _loop1min, int _loop1max, int _loop2min, int _loop2max, double _distanceStem1, double _distanceStem2, double _stemRmsdTol, double _totalRmsdTol, bool _matchFirstStemOnly){

        int numDualLoops = 0;

        if (_stemResidues1.size() != 2) throw MslSizeException(MslTools::stringf("PDBFragments::searchForMatchingDualFragments stemResidue1 != 2, it = %d",_stemResidues1.size()));

        Position &stem1res1 = _sys1.getPosition(_stemResidues1[0]);
        Position &stem1res2 = _sys1.getPosition(_stemResidues1[1]);

	if (! (stem1res1.atomExists("CA") && stem1res1.atomExists("N") && stem1res1.atomExists("C"))) {
	  throw MslAtomNotExistException(MslTools::stringf("PDBFragments::searchForMatchingDualFragments one of N-CA-C doesn't exist on %s",stem1res1.getPositionId().c_str()));
	};

	if (! (stem1res2.atomExists("CA") && stem1res2.atomExists("N") && stem1res2.atomExists("C"))) {
	  throw MslAtomNotExistException(MslTools::stringf("PDBFragments::searchForMatchingDualFragments one of N-CA-C doesn't exist on %s",stem1res2.getPositionId().c_str()));
	};

	AtomContainer stem1;
	stem1.addAtom(stem1res1.getAtom("CA"));
	stem1.addAtom(stem1res2.getAtom("CA"));

	
        if (_stemResidues2.size() != 2) throw MslSizeException(MslTools::stringf("PDBFragments::searchForMatchingDualFragments stemResidue2 != 2, it = %d",_stemResidues2.size()));

        Position &stem2res1 = _sys2.getPosition(_stemResidues2[0]);
        Position &stem2res2 = _sys2.getPosition(_stemResidues2[1]);

	if (! (stem2res1.atomExists("CA") && stem2res1.atomExists("N") && stem2res1.atomExists("C"))) {
	  throw MslAtomNotExistException(MslTools::stringf("PDBFragments::searchForMatchingDualFragments one of N-CA-C doesn't exist on %s",stem2res1.getPositionId().c_str()));
	};
	if (! (stem2res2.atomExists("CA") && stem2res2.atomExists("N") && stem2res2.atomExists("C"))) {
	  throw MslAtomNotExistException(MslTools::stringf("PDBFragments::searchForMatchingDualFragments one of N-CA-C doesn't exist on %s",stem2res2.getPositionId().c_str()));
	};	

	AtomContainer stem2;
	stem2.addAtom(stem2res1.getAtom("CA"));
	stem2.addAtom(stem2res2.getAtom("CA"));

	AtomContainer stemsFullBB;
	stemsFullBB.addAtom(stem1res1.getAtom("N"));
	stemsFullBB.addAtom(stem1res1.getAtom("CA"));
	stemsFullBB.addAtom(stem1res1.getAtom("C"));
	stemsFullBB.addAtom(stem2res1.getAtom("N"));
	stemsFullBB.addAtom(stem2res1.getAtom("CA"));
	stemsFullBB.addAtom(stem2res1.getAtom("C"));	
	stemsFullBB.addAtom(stem1res2.getAtom("N"));
	stemsFullBB.addAtom(stem1res2.getAtom("CA"));
	stemsFullBB.addAtom(stem1res2.getAtom("C"));
	stemsFullBB.addAtom(stem2res2.getAtom("N"));
	stemsFullBB.addAtom(stem2res2.getAtom("CA"));
	stemsFullBB.addAtom(stem2res2.getAtom("C"));	

	// Define distances between the stems on each structure
	if (_distanceStem1 == 0.0){
	  _distanceStem1 = stem1res1.getAtom("CA").distance2(stem1res2.getAtom("CA")); 
	}	
	if (_distanceStem2 == 0.0){
	  _distanceStem2 = stem2res1.getAtom("CA").distance2(stem2res2.getAtom("CA")); 
	}	

	// Fragment1: r1 to r2
	// Fragment2: r3 to r4

	// Each chain
	for (uint r1 = 0 ; r1 < fragDB.size()-_loop1min+_loop2min;r1++){	

	       Atom &loop1Res1  = fragDB(r1);

	       for (uint r2 = r1+_loop1min; r2 < r1+_loop1max && r2 < fragDB.size()-_loop1min+_loop2min;r2++){
		 
		 Atom &loop1Res2  = fragDB(r2);

		 // Same PDB, Same Chain, please..
		 if (loop1Res1.getSegID() != loop1Res2.getSegID()) break;
		 if (loop1Res1.getChainId() != loop1Res2.getChainId()) break;

		 // Look for second fragment..
		 for (uint r3 = r2+1; r3 < fragDB.size()-_loop2min;r3++){

		   Atom &loop2Res1  = fragDB(r3);

		   // Same PDB, Same Chain, please..
		   if (loop1Res1.getSegID()   != loop2Res1.getSegID()) break;
		   if (loop1Res1.getChainId() != loop2Res1.getChainId()) break;

		   // Look for R2-R3 distance match
		   double dist2 = loop1Res2.distance2(loop2Res1);
		   if (fabs(dist2 - _distanceStem2) > 4) {
		     continue;
		   }

		   for (uint r4 = r3+_loop2min; r4 < r3+_loop2max && r4 < fragDB.size();r4++){
		     if (r1 == r4 || r1 == r3) continue;

		     Atom &loop2Res2  = fragDB(r4);

		     // Same PDB, Same Chain, please..
		     if (loop1Res1.getSegID()   != loop2Res2.getSegID()) break;
		     if (loop1Res1.getChainId() != loop2Res2.getChainId()) break;

		     // Look for R1-R4 distance match (< 1 Angstroms)
		     double dist1 = loop1Res1.distance2(loop2Res2);
		     if (fabs(dist1 - _distanceStem1) > 4){
		       continue;
		     } 

		     // RMSD Check to stems
		     AtomContainer testMatchStem1;
		     testMatchStem1.addAtom(loop1Res1);
		     testMatchStem1.addAtom(loop2Res2);

		     AtomContainer testMatchStem2;
		     testMatchStem2.addAtom(loop1Res2);
		     testMatchStem2.addAtom(loop2Res1);
		       

		     AtomPointerVector matchStems;
		     AtomPointerVector refStems;

		     Transforms t;
		     if (_matchFirstStemOnly){
		       // Match stems independently, the "matchStems" only has stem1 in it.

		       // Move testMatchStem1+2 onto stem1 , using testMatchStem1 onto stem1.
		       testMatchStem1.saveCoor("pre");
		       if (!t.rmsdAlignment(testMatchStem1.getAtomPointers(),stem1.getAtomPointers()))
			 throw MslRmsdAlignmentFailException("PDBFragments::searchForMatchingDualFragments testMatchStem1,stem1");

		       double rmsd = testMatchStem1.getAtomPointers().rmsd(stem1.getAtomPointers());

		       if (rmsd > _stemRmsdTol){
			 MSLOUT.stream() << "RMSD TOO LARGE: "<<rmsd<<" tol: "<<_stemRmsdTol<<" on stem1."<<endl;
			 continue;
		       } 


		       testMatchStem2.saveCoor("pre");
		       if (!t.rmsdAlignment(testMatchStem2.getAtomPointers(),stem2.getAtomPointers()))
			 throw MslRmsdAlignmentFailException("PDBFragments::searchForMatchingDualFragments testMatchStem2,stem2");

		       rmsd = testMatchStem2.getAtomPointers().rmsd(stem2.getAtomPointers());
		       if (rmsd > _stemRmsdTol){
			 MSLOUT.stream() << "RMSD TOO LARGE: "<<rmsd<<" tol: "<<_stemRmsdTol<<" on stem2."<<endl;
			 continue;
		       } 		       
		       matchStems = testMatchStem1.getAtomPointers();
		       refStems   = stem1.getAtomPointers();

		     } else {

		       matchStems = testMatchStem1.getAtomPointers() + testMatchStem2.getAtomPointers();
		       refStems = stem1.getAtomPointers() + stem2.getAtomPointers();

		       // Move testMatchStem1+2 onto stem1 , using testMatchStem1 onto stem1.
		       matchStems.saveCoor("pre");
		       if (!t.rmsdAlignment(matchStems,refStems)) 
			 throw MslRmsdAlignmentFailException("PDBFragments::searchForMatchingDualFragments matchStems, refStems, dualLoops)");

		       double rmsd = matchStems.rmsd(refStems);
		       if (rmsd > _stemRmsdTol){
			 //MSLOUT.stream() << "RMSD TOO LARGE: "<<rmsd<<" tol: "<<_stemRmsdTol<<" on dualLoop stems.("<<loop1Res1.getSegID()<<")"<<endl;
			 continue;
		       }

		     }

		     // IF nothing further, simply printout the match
		     if (pdbDir == ""){
		       MSLOUT.fprintf(stdout,"MATCH %s %s-%s %s-%s\n",loop1Res1.getSegID().c_str(), loop1Res1.getPositionId().c_str(),loop1Res2.getPositionId().c_str(),loop2Res1.getPositionId().c_str(),loop2Res2.getPositionId().c_str());
		       continue;
		     }


		     // Load PDB and extract region
		     string allAtomFileName = MslTools::stringf("%s/%s.pdb",pdbDir.c_str(),loop1Res1.getSegID().c_str());
		     
		     System allAtomSys;
		     allAtomSys.readPdb(allAtomFileName);

		     // Reset segid..
		     for (uint ats = 0; ats < allAtomSys.getAtomPointers().size();ats++){
		       allAtomSys.getAtom(ats).setSegID("");
		     }

		     // Good stem1,stem2 alignment, move the loops onto 
		     AtomContainer *dualLoops = new AtomContainer();

		     for (uint m = r1; m <= r2;m++){
		       (*dualLoops).addAtoms(allAtomSys.getPosition(fragDB(m).getPositionId()).getAtomPointers());
		     }

		     for (uint m = r3; m <= r4;m++){
		       (*dualLoops).addAtoms(allAtomSys.getPosition(fragDB(m).getPositionId()).getAtomPointers());
		     }


		     AtomContainer fullBB;
		     fullBB.addAtom(allAtomSys.getPosition(fragDB(r1).getPositionId()).getAtom("N"));
		     fullBB.addAtom(allAtomSys.getPosition(fragDB(r1).getPositionId()).getAtom("CA"));
		     fullBB.addAtom(allAtomSys.getPosition(fragDB(r1).getPositionId()).getAtom("C"));
		     fullBB.addAtom(allAtomSys.getPosition(fragDB(r2).getPositionId()).getAtom("N"));
		     fullBB.addAtom(allAtomSys.getPosition(fragDB(r2).getPositionId()).getAtom("CA"));
		     fullBB.addAtom(allAtomSys.getPosition(fragDB(r2).getPositionId()).getAtom("C"));
		     fullBB.addAtom(allAtomSys.getPosition(fragDB(r3).getPositionId()).getAtom("N"));
		     fullBB.addAtom(allAtomSys.getPosition(fragDB(r3).getPositionId()).getAtom("CA"));
		     fullBB.addAtom(allAtomSys.getPosition(fragDB(r3).getPositionId()).getAtom("C"));
		     fullBB.addAtom(allAtomSys.getPosition(fragDB(r4).getPositionId()).getAtom("N"));
		     fullBB.addAtom(allAtomSys.getPosition(fragDB(r4).getPositionId()).getAtom("CA"));
		     fullBB.addAtom(allAtomSys.getPosition(fragDB(r4).getPositionId()).getAtom("C"));		     

		     // Move the dualLoops onto the correct frame		       
		     if (!t.rmsdAlignment(fullBB.getAtomPointers(),stemsFullBB.getAtomPointers(),(*dualLoops).getAtomPointers())) 
			 throw MslRmsdAlignmentFailException("PDBFragments::searchForMatchingDualFragments matchStems, refStems, dualLoops)");
		     if (!t.rmsdAlignment(fullBB.getAtomPointers(),stemsFullBB.getAtomPointers())) 
			 throw MslRmsdAlignmentFailException("PDBFragments::searchForMatchingDualFragments matchStems, refStems, dualLoops)");

		     double rmsd = fullBB.getAtomPointers().rmsd(stemsFullBB.getAtomPointers());
		     if (rmsd > _totalRmsdTol) {
		       MSLOUT.stream() << "BB-RMSD too large: "<<rmsd<<" tol: "<<_totalRmsdTol<<endl;
		       continue;
		     }
		     
		     lastResults.push_back(dualLoops);		     
		     
		     MSLOUT.fprintf(stdout,"MATCH %5d %s %s-%s %s-%s\n",numDualLoops,loop1Res1.getSegID().c_str(), loop1Res1.getPositionId().c_str(),loop1Res2.getPositionId().c_str(),loop2Res1.getPositionId().c_str(),loop2Res2.getPositionId().c_str());
		     numDualLoops++;

		   } // FOR r4
		 } // FOR r3
	       } // FOR r2
	     } // FOR r1

	return numDualLoops;

}
int PDBFragments::searchForMatchingDualFragments(System &_sys1, std::vector<std::string> &_stemResidues1,
						 System &_sys2, std::vector<std::string> &_stemResidues2,
						 System &_searchSys, int _loop1min, int _loop1max, int _loop2min, int _loop2max, double _distanceStem1, double _distanceStem2, double _stemRmsdTol, double _totalRmsdTol, bool _matchFirstStemOnly){

        int numDualLoops = 0;

        if (_stemResidues1.size() != 2) throw MslSizeException(MslTools::stringf("PDBFragments::searchForMatchingDualFragments stemResidue1 != 2, it = %d",_stemResidues1.size()));

        Position &stem1res1 = _sys1.getPosition(_stemResidues1[0]);
        Position &stem1res2 = _sys1.getPosition(_stemResidues1[1]);

	if (! (stem1res1.atomExists("CA") && stem1res1.atomExists("N") && stem1res1.atomExists("C"))) {
	  throw MslAtomNotExistException(MslTools::stringf("PDBFragments::searchForMatchingDualFragments one of N-CA-C doesn't exist on %s",stem1res1.getPositionId().c_str()));
	};

	if (! (stem1res2.atomExists("CA") && stem1res2.atomExists("N") && stem1res2.atomExists("C"))) {
	  throw MslAtomNotExistException(MslTools::stringf("PDBFragments::searchForMatchingDualFragments one of N-CA-C doesn't exist on %s",stem1res2.getPositionId().c_str()));
	};

	AtomContainer stem1;
	stem1.addAtom(stem1res1.getAtom("N"));
	stem1.addAtom(stem1res1.getAtom("CA"));
	stem1.addAtom(stem1res1.getAtom("C"));
	stem1.addAtom(stem1res2.getAtom("N"));
	stem1.addAtom(stem1res2.getAtom("CA"));
	stem1.addAtom(stem1res2.getAtom("C"));	   

        if (_stemResidues2.size() != 2) throw MslSizeException(MslTools::stringf("PDBFragments::searchForMatchingDualFragments stemResidue2 != 2, it = %d",_stemResidues2.size()));

        Position &stem2res1 = _sys2.getPosition(_stemResidues2[0]);
        Position &stem2res2 = _sys2.getPosition(_stemResidues2[1]);

	
	if (! (stem2res1.atomExists("CA") && stem2res1.atomExists("N") && stem2res1.atomExists("C"))) {
	  throw MslAtomNotExistException(MslTools::stringf("PDBFragments::searchForMatchingDualFragments one of N-CA-C doesn't exist on %s",stem2res1.getPositionId().c_str()));
	};
	if (! (stem2res2.atomExists("CA") && stem2res2.atomExists("N") && stem2res2.atomExists("C"))) {
	  throw MslAtomNotExistException(MslTools::stringf("PDBFragments::searchForMatchingDualFragments one of N-CA-C doesn't exist on %s",stem2res2.getPositionId().c_str()));
	};	

	AtomContainer stem2;
	stem2.addAtom(stem2res1.getAtom("N"));
	stem2.addAtom(stem2res1.getAtom("CA"));
	stem2.addAtom(stem2res1.getAtom("C"));
	stem2.addAtom(stem2res2.getAtom("N"));
	stem2.addAtom(stem2res2.getAtom("CA"));
	stem2.addAtom(stem2res2.getAtom("C"));	   

	AtomPointerVector stems = stem1.getAtomPointers() + stem2.getAtomPointers();
	// Define distances between the stems on each structure
	if (_distanceStem1 == 0.0){
	  _distanceStem1 = stem1res1.getAtom("CA").distance2(stem1res2.getAtom("CA")); 
	}	
	if (_distanceStem2 == 0.0){
	  _distanceStem2 = stem2res1.getAtom("CA").distance2(stem2res2.getAtom("CA")); 
	}	

	// Fragment1: r1 to r2
	// Fragment2: r3 to r4

	// Each chain
	for (uint c = 0; c < _searchSys.chainSize();c++){
	     
	     // Walk through residues
	     Chain &ch = _searchSys.getChain(c);

	     for (uint r1 = 0; r1 < ch.positionSize();r1++){
	       Residue &loop1Res1  = ch.getResidue(r1);

	       // Skip over residues with out a proper backbone
	       if (! (loop1Res1.atomExists("CA") && loop1Res1.atomExists("N") && loop1Res1.atomExists("C"))) continue;


	       for (uint r2 = r1+_loop1min; r2 < r1+_loop1max && r2 < ch.positionSize();r2++){
		 Residue &loop1Res2  = ch.getResidue(r2);

		 // Skip over residues without a proper backbone
		 if (! (loop1Res2.atomExists("CA") && loop1Res2.atomExists("N") && loop1Res2.atomExists("C"))) continue;

		 // Look for second fragment..
		 for (uint r3 = r2+1; r3 < ch.positionSize();r3++){
		   if (r2 == r3) continue;

		   Residue &loop2Res1  = ch.getResidue(r3);

		   // Skip pover residues without a proper backbone
		   if (! (loop2Res1.atomExists("CA") && loop2Res1.atomExists("N") && loop2Res1.atomExists("C"))) continue;

		   // Look for R2-R3 distance match
		   double dist2 = loop1Res2.distance(loop2Res1,"CA",true); // true = distance squared
		   if (fabs(dist2 - _distanceStem2) > 20) {
		     continue;
		   }

		   for (uint r4 = r3+_loop2min; r4 < r3+_loop2max && r4 < ch.positionSize();r4++){
		     if (r1 == r4 || r1 == r3) continue;

		     Residue &loop2Res2  = ch.getResidue(r4);
		     if (! (loop2Res2.atomExists("CA") && loop2Res2.atomExists("N") && loop2Res2.atomExists("C"))) continue;

		     // Look for R1-R4 distance match (< 1 Angstroms)
		     double dist1 = loop1Res1.distance(loop2Res2,"CA",true); // true = distance squared
		     if (fabs(dist1 - _distanceStem1) > 20){
		       continue;
		     } 


		     // RMSD Check to stems
		     AtomContainer testMatchStem1;
		     testMatchStem1.addAtom(loop1Res1("N"));
		     testMatchStem1.addAtom(loop1Res1("CA"));
		     testMatchStem1.addAtom(loop1Res1("C"));
		     testMatchStem1.addAtom(loop2Res2("N"));
		     testMatchStem1.addAtom(loop2Res2("CA"));
		     testMatchStem1.addAtom(loop2Res2("C"));


		     Transforms t;

		     // Move testMatchStem1+2 onto stem1 , using testMatchStem1 onto stem1.
		     testMatchStem1.saveCoor("pre");
		     if (!t.rmsdAlignment(testMatchStem1.getAtomPointers(),stem1.getAtomPointers()))
		       throw MslRmsdAlignmentFailException("PDBFragments::searchForMatchingDualFragments testMatchStem1,stem1");

		     double rmsd = testMatchStem1.getAtomPointers().rmsd(stem1.getAtomPointers());

		     if (rmsd > _stemRmsdTol){
		       MSLOUT.stream() << "RMSD TOO LARGE: "<<rmsd<<" tol: "<<_stemRmsdTol<<" on stem1."<<endl;
		       continue;
		     } 

		     AtomContainer testMatchStem2;
		     testMatchStem2.addAtom(loop1Res2("N"));
		     testMatchStem2.addAtom(loop1Res2("CA"));
		     testMatchStem2.addAtom(loop1Res2("C"));
		     testMatchStem2.addAtom(loop2Res1("N"));
		     testMatchStem2.addAtom(loop2Res1("CA"));
		     testMatchStem2.addAtom(loop2Res1("C"));

		     testMatchStem2.saveCoor("pre");
		     if (!t.rmsdAlignment(testMatchStem2.getAtomPointers(),stem2.getAtomPointers()))
		       throw MslRmsdAlignmentFailException("PDBFragments::searchForMatchingDualFragments testMatchStem2,stem2");

		     rmsd = testMatchStem2.getAtomPointers().rmsd(stem2.getAtomPointers());
		     if (rmsd > _stemRmsdTol){
		       MSLOUT.stream() << "RMSD TOO LARGE: "<<rmsd<<" tol: "<<_stemRmsdTol<<" on stem2."<<endl;
		       continue;
		     } 
		       
		     // Good stem1,stem2 alignment, move the loops onto 
		     AtomContainer *dualLoops = new AtomContainer();
		     for (uint m = r1; m <= r2;m++){
		       (*dualLoops).addAtoms(ch.getResidue(m).getAtomPointers());
		     }

		     for (uint m = r3; m <= r4;m++){
		       (*dualLoops).addAtoms(ch.getResidue(m).getAtomPointers());
		     }


		     // Move the dualLoops onto the correct frame		       
		     testMatchStem1.applySavedCoor("pre");
		     AtomPointerVector matchStems;
		     AtomPointerVector refStems;
		     if (_matchFirstStemOnly){
		       matchStems = testMatchStem1.getAtomPointers();
		       refStems = stem1.getAtomPointers();
		     } else {
		       testMatchStem2.applySavedCoor("pre");
		       matchStems = testMatchStem1.getAtomPointers() + testMatchStem2.getAtomPointers();
		       refStems = stem1.getAtomPointers() + stem2.getAtomPointers();
		     }

		     if (!t.rmsdAlignment(matchStems,refStems,(*dualLoops).getAtomPointers())) 
			 throw MslRmsdAlignmentFailException("PDBFragments::searchForMatchingDualFragments matchStems, refStems, dualLoops)");

		     rmsd = matchStems.rmsd(refStems);
		     if (rmsd > _totalRmsdTol){
		       MSLOUT.stream() << "RMSD TOO LARGE: "<<rmsd<<" tol: "<<_totalRmsdTol<<" on dualLoop stems."<<endl;
		       delete(dualLoops);
		       continue;
		     }

		     lastResults.push_back(dualLoops);		     
		     
		     
		     numDualLoops++;

		   } // FOR r4
		 } // FOR r3
	       } // FOR r2
	     } // FOR r1
	} // For chain c

	return numDualLoops;
}

int PDBFragments::searchForMatchingFragmentsSpots(System &_sys, std::vector<std::string> &_stemResidues, int _maxResiduesBetweenStems, double _rmsdTol){
	if (fragDB.size() <= 4){
		cerr << "ERROR 1232 PDBFragments::searchForMatchingFragments(), fragment database is too small, maybe you forgot to load it using loadFragmentDatabase() ?"<<endl;
		exit(1232);
	}
	int numFrags = 0;

	// Remove last set of results
	lastResults.clear();

	// Clear our set of matched sequences
	matchedSequences.clear();


	if (_stemResidues.size() != 3){
	  cerr << "ERROR 249582 PDBFragments::searchForMatchingFragments3, need exactly 3 stem residues, got: "<<_stemResidues.size()<<endl;
	  exit(249582);
	}

	if (!_sys.positionExists(_stemResidues[0]) ||!_sys.positionExists(_stemResidues[1]) ||!_sys.positionExists(_stemResidues[2])){
	  cerr << "ERROR 23423 PDBFragments::searchForMatchingFragments3, couldn't find a stem residue in pdb file: "<<_stemResidues[0]<<" , "<<_stemResidues[1]<<" , "<<_stemResidues[2]<<endl;
	  exit(23423);
	}
	Position &pos1 = _sys.getPosition(_stemResidues[0]);
	Position &pos2 = _sys.getPosition(_stemResidues[1]);
	Position &pos3 = _sys.getPosition(_stemResidues[2]);

	AtomPointerVector stemBBats;
	stemBBats.push_back(&pos1.getAtom("N"));
	stemBBats.push_back(&pos1.getAtom("CA"));
	stemBBats.push_back(&pos1.getAtom("C"));

	stemBBats.push_back(&pos2.getAtom("N"));
	stemBBats.push_back(&pos2.getAtom("CA"));
	stemBBats.push_back(&pos2.getAtom("C"));

	stemBBats.push_back(&pos3.getAtom("N"));
	stemBBats.push_back(&pos3.getAtom("CA"));
	stemBBats.push_back(&pos3.getAtom("C"));

         
        // Store stem-to-stem distance-squared vector (for filtering candidates below)
        AtomPointerVector stems;
	stems.push_back(&pos1.getAtom("CA"));
	stems.push_back(&pos2.getAtom("CA"));
	stems.push_back(&pos3.getAtom("CA"));
	
	vector<double> stemDistanceSq;
	for (uint c = 0; c < stems.size();c++){
	  for (uint n = c+1; n < stems.size();n++){
	    double distSq = stems(c).distance2(stems(n));			
	    stemDistanceSq.push_back(distSq);
	  }
	}

	MSLOUT.stream() << "FragDB.size(): "<<fragDB.size()<<endl;


	// Now loop over all fragments in database checking for ones of the correct size
	double tol = 64; // Tolerance of distance to be deviant from stems in Angstroms^2
	int matchIndex = 0; // Index for keeping track of matches
	for (uint i = 0 ; i < fragDB.size()-(2*_maxResiduesBetweenStems+stems.size());i++){


			//bool validTriplet = true;
			int secondPositionIndex = -1;

			// Scan ahead for second position
			int closestDistanceIndex = 0;
			double minDiff           = MslTools::doubleMax;
			for (uint j = 1; j <= _maxResiduesBetweenStems;j++){

			  // Second position is no good due to segid match
			  if (fragDB[i]->getSegID() != fragDB[i+j]->getSegID()){
			    //validTriplet = false;
			    break;
			  }

			  // Second position is no good due to chain match
			  if (fragDB[i]->getChainId() != fragDB[i+j]->getChainId()){
			    //validTriplet = false;
			    break;
			  }

			  // Second position is no good due to gap..
			  if ( abs(fragDB[i]->getResidueNumber() - fragDB[i+j]->getResidueNumber()) > j){
			    //validTriplet = false;
			    break;
			  }


			  // Second position is no good due to distance tolerance.., don't break loop!
			  double distSq = fragDB[i]->distance2(*fragDB[i+j]);
			  if (fabs(distSq - stemDistanceSq[0]) < minDiff){
			    minDiff              = fabs(distSq - stemDistanceSq[0]);
			    closestDistanceIndex = j;
			  }
			}

			if (minDiff > tol) continue;
			secondPositionIndex = i + closestDistanceIndex;

			if (secondPositionIndex+_maxResiduesBetweenStems > fragDB.size()) continue;
			int thirdPositionIndex = -1;
			closestDistanceIndex = 0;
			minDiff = MslTools::doubleMax;
			double closestDistance2   = MslTools::doubleMax;
			double closestDistance3   = MslTools::doubleMax;
			for (uint j = 1; j <= _maxResiduesBetweenStems;j++){

			  // Second position is no good due to segid match
			  if (fragDB[secondPositionIndex]->getSegID() != fragDB[secondPositionIndex+j]->getSegID()){
			    //validTriplet = false;
			    break;
			  }

			  // Second position is no good due to chain match
			  if (fragDB[secondPositionIndex]->getChainId() != fragDB[secondPositionIndex+j]->getChainId()){
			    //validTriplet = false;
			    break;
			  }

			  // Second position is no good due to gap..
			  if ( abs(fragDB[secondPositionIndex]->getResidueNumber() - fragDB[secondPositionIndex+j]->getResidueNumber()) > j){
			    //validTriplet = false;
			    break;
			  }


			  // Second position is no good due to distance tolerance.., don't break loop!
			  double distSq2 = fragDB[i]->distance2(*fragDB[secondPositionIndex+j]);
			  double distSq3 = fragDB[secondPositionIndex]->distance2(*fragDB[secondPositionIndex+j]);
			  double diffFound = fabs(distSq2 - stemDistanceSq[1]) + fabs(distSq3 - stemDistanceSq[2]);
			  if (diffFound < minDiff){
			    closestDistanceIndex = j;
			    closestDistance2      = distSq2;
			    closestDistance3      = distSq3;
			    minDiff               = diffFound;
			  }

			}			  

			if (fabs(closestDistance2-stemDistanceSq[1]) > tol && fabs(closestDistance3 - stemDistanceSq[2]) > tol) continue;
			thirdPositionIndex = secondPositionIndex + closestDistanceIndex;

			if (thirdPositionIndex >= fragDB.size()) { cerr << "ERROR 3rd position index too large"<<endl; continue;}
			
			// MSLOUT.stream() << MslTools::stringf("%06d %s-%1s_%04d%s ",
			// 				     matchIndex,
			// 				     fragDB(i).getSegID().c_str(),
			// 				     fragDB(i).getChainId().c_str(),
			// 				     fragDB(i).getResidueNumber(),
			// 				     fragDB(i).getResidueIcode().c_str());
			// MSLOUT.stream() << MslTools::stringf(" %s-%1s_%04d%s ",
			// 			       fragDB(secondPositionIndex).getSegID().c_str(),
			// 			       fragDB(secondPositionIndex).getChainId().c_str(),
			// 			       fragDB(secondPositionIndex).getResidueNumber(),
  			// 			       fragDB(secondPositionIndex).getResidueIcode().c_str());
			// MSLOUT.stream() << MslTools::stringf(" %s-%1s_%04d%s \n",
			// 			       fragDB(thirdPositionIndex).getSegID().c_str(),
			// 			       fragDB(thirdPositionIndex).getChainId().c_str(),
			// 			       fragDB(thirdPositionIndex).getResidueNumber(),
			// 			       fragDB(thirdPositionIndex).getResidueIcode().c_str());


			// align and print winning fragment
			// Make a copy of frag stem, so not to effect original atoms
			AtomPointerVector fragStem; 
			fragStem.push_back(fragDB[i]);
			fragStem.push_back(fragDB[secondPositionIndex]);
			fragStem.push_back(fragDB[thirdPositionIndex]);
			// Align the fragStem to stem to see if it matches well enough...
			Transforms tm;
			fragStem.saveCoor("pre");
			tm.rmsdAlignment(fragStem,stems);

			double rmsd = fragStem.rmsd(stems);


			// Continue if the RMSD filter not passed
			if (rmsd > 0.5){
				continue;
			}

			matchIndex++;
			fragStem.applySavedCoor("pre");

			// Create new atoms for middle atoms..
			AtomContainer fullFragCA;
			string matchSeq = "";
			for (uint f = i; f <= thirdPositionIndex;f++){
			        fullFragCA.addAtom(*fragDB[f]);
				matchSeq += MslTools::getOneLetterCode(fragDB[f]->getResidueName());
			}



			string key = MslTools::stringf("%06d-%s-%1s_%04d%s-%1s_%04d%s",
						       matchIndex,
						       fragDB(i).getSegID().c_str(),
						       fragDB(i).getChainId().c_str(),
						       fragDB(i).getResidueNumber(),
						       fragDB(i).getResidueIcode().c_str(),
						       fragDB(thirdPositionIndex).getChainId().c_str(),
						       fragDB(thirdPositionIndex).getResidueNumber(),
						       fragDB(thirdPositionIndex).getResidueIcode().c_str());

			matchedSequences[key] = matchSeq;

			// if (_regex != ""){
			//   if (!boost::regex_search(matchSeq.c_str(),boost::regex(_regex))){
			//     continue;
			//   } else {
			//     MSLOUT.stream() << "RegEx Matched."<<endl;
			//   }
			// }
			

			
			tm.rmsdAlignment(fragStem,stems,fullFragCA.getAtomPointers());




			
			bool successful = true;

			  if (pdbDir != ""){

			              Atom &at1 = fullFragCA.getAtom(0);
			              Atom &at2 = fullFragCA.getAtom(fullFragCA.size()-1);

				      string allAtomFileName = MslTools::stringf("%s/%s.pdb",pdbDir.c_str(),at1.getSegID().c_str());

				      System allAtomSys;
				      allAtomSys.readPdb(allAtomFileName);
				      for (uint ats = 0; ats < allAtomSys.getAtomPointers().size();ats++){
					allAtomSys.getAtom(ats).setSegID("");
				      }

				      AtomSelection sel2(allAtomSys.getAtomPointers());

			    	      stringstream ss;
				      char tmpstr[100];
				      //sprintf(tmpstr,"chain %1s and resi %d-%-d and name CA",at1.getChainId().c_str(),at1.getResidueNumber(),at2.getResidueNumber());
				      sprintf(tmpstr,"chain %1s and resi %d-%-d and name CA",	fragDB[i]->getChainId().c_str(),fragDB[i]->getResidueNumber(),fragDB[thirdPositionIndex]->getResidueNumber());
				      ss << tmpstr;
				      AtomPointerVector caAts = sel2.select(ss.str());


				      if (!tm.rmsdAlignment(caAts,fullFragCA.getAtomPointers(),allAtomSys.getAtomPointers())){
					MSLOUT.stream() << "Problem aligning all atoms using the C-alpha trace"<<endl;
					MSLOUT.stream() << "PDB: "<<allAtomFileName<<endl;
					MSLOUT.stream() << "\tTrying to align with: "<<ss.str()<<endl;
					MSLOUT.stream() << "\tSystem: "<<allAtomSys.getSizes();
					MSLOUT.stream() << "\tSelected: "<<caAts.size()<<" atoms using '"<<ss.str()<<"'"<<endl;
					MSLOUT.stream() << "\tReference: "<<fullFragCA.getAtomPointers().size()<<" atoms"<<endl;
					MSLOUT.stream() << fullFragCA.getAtomPointers();
					continue;
				      } 

				      // Get BB atoms from stem-equivalent residues
				      AtomPointerVector allAts_stemBBats;
				      Position &stem1eq  = allAtomSys.getPosition(MslTools::getPositionId(fragDB[i]->getChainId(),fragDB[i]->getResidueNumber(),fragDB[i]->getResidueIcode()));
				      allAts_stemBBats.push_back(&stem1eq.getAtom("N"));
				      allAts_stemBBats.push_back(&stem1eq.getAtom("CA"));
				      allAts_stemBBats.push_back(&stem1eq.getAtom("C"));

				      Position &stem2eq  = allAtomSys.getPosition(MslTools::getPositionId(fragDB[secondPositionIndex]->getChainId(),fragDB[secondPositionIndex]->getResidueNumber(),fragDB[secondPositionIndex]->getResidueIcode()));
				      allAts_stemBBats.push_back(&stem2eq.getAtom("N"));
				      allAts_stemBBats.push_back(&stem2eq.getAtom("CA"));
				      allAts_stemBBats.push_back(&stem2eq.getAtom("C"));

				      Position &stem3eq  = allAtomSys.getPosition(MslTools::getPositionId(fragDB[thirdPositionIndex]->getChainId(),fragDB[thirdPositionIndex]->getResidueNumber(),fragDB[thirdPositionIndex]->getResidueIcode()));
				      allAts_stemBBats.push_back(&stem3eq.getAtom("N"));
				      allAts_stemBBats.push_back(&stem3eq.getAtom("CA"));
				      allAts_stemBBats.push_back(&stem3eq.getAtom("C"));
				      	
				      double bbRMSD = allAts_stemBBats.rmsd(stemBBats);


				      if (bbRMSD > _rmsdTol){
					continue;
				      }
				      fprintf(stdout,"(%4s and chain %1s and resi %3d-%3d)  %8.3f ",
					      fragDB[i]->getSegID().c_str(),
					      fragDB[i]->getChainId().c_str(),
					      fragDB[i]->getResidueNumber(),
					      fragDB[thirdPositionIndex]->getResidueNumber(),
					      rmsd);

				      fprintf(stdout, "%8.3f",bbRMSD);

				      ss.str("");
				      char tmpstr2[100];
				      sprintf(tmpstr2,"chain %1s and resi %d-%-d",at1.getChainId().c_str(),at1.getResidueNumber(),at2.getResidueNumber());
				      ss << tmpstr2;
				      AtomPointerVector allAts = sel2.select(ss.str());
				      lastResults.push_back(new AtomContainer());
				      lastResults.back()->addAtoms(allAts);




			  }
	 
			fprintf(stdout,"\n");
			fragStem.applySavedCoor("pre");

			if (successful){
				numFrags++;
			}
			
	}
	MSLOUT.stream() << "Done."<<endl;
	
	return numFrags;
}


int PDBFragments::searchForMatchingFragmentsLinear(System &_sys, string &_startRes, string &_endRes, string _regex, double _rmsdTol, int _maxFrags){

	if (  ! ( _sys.positionExists(_startRes)  && _sys.positionExists(_endRes) ) ) {
		cerr << "ERROR 2342 residue(s) don't exist: "<<_startRes<<" "<<_endRes<<endl;
		return -1;
	}

	int numFrags = 0;
  
	// Remove last set of results
	lastResults.clear();

	// Clear our set of matched sequences
	matchedSequences.clear();

	// AtomVector of CA/backbone atoms
	int residueSeparation = _sys.getPositionIndex(_endRes) - _sys.getPositionIndex(_startRes);
	AtomPointerVector bbAts;
	for (uint i = _sys.getPositionIndex(_startRes); i <= _sys.getPositionIndex(_endRes); i++){
		if (_sys.getPosition(i).atomExists("CA")){
			bbAts.push_back(&_sys.getPosition(i).getAtom("CA"));
		} else {
			cerr << "Position: "<<_sys.getPosition(i).toString()<<" doesn't have a CA atom"<<endl;
			return -1;
		}
	}

	MSLOUT.stream() << "FragDB.size(): "<<fragDB.size()<< " "<<residueSeparation<<" ."<<_regex<<"."<<endl;
	// Now loop over all fragments in database checking for ones of the correct size
	int matchIndex = 0; // Index for keeping track of matches
	Transforms tm;
	stringstream ss;
	for (uint i = 0 ; i < fragDB.size()-residueSeparation; i++){
		//MSLOUT.stream() << "Fragment("<<i<<") "<<endl;

		// Break if we already found enough fragments
		if (_maxFrags != -1 && numFrags > _maxFrags){
			break;
		}
	
		// Filter by pdb/chain breaks
		if (fragDB(i).getSegID() != fragDB(i+residueSeparation).getSegID()){
			//MSLOUT.stream() << "PDB break: "<<fragDB(i).getSegID() << " "<<fragDB(i+residueSeparation).getSegID()<<endl;
			continue;
		}
		if (fragDB(i).getChainId() != fragDB(i+residueSeparation).getChainId()){
			//MSLOUT.stream() << "CHAIN break"<<endl;
			continue;
		}

		AtomContainer fragBB;
		string matchSeq = "";
		for (uint j = 0; j <= residueSeparation; j++){
			fragBB.addAtom(fragDB(i+j));
			matchSeq += MslTools::getOneLetterCode(fragDB(i+j).getResidueName());
		}

		//fragBB.saveCoor("pre");
		if (!tm.rmsdAlignment(fragBB.getAtomPointers(), bbAts)){
			cerr << "ERROR in alignment: "<<fragBB.size()<<" "<<bbAts.size()<<endl;
			continue;
		}
		// Simple distance squared check to speed things up.
		//MSLOUT.stream() << "Distance fragBB(0) to bbAts(0): " << fragBB(0).distance(bbAts(0))<<" fragBB:"<<fragBB(0).toString()<<" "<<fragBB(fragBB.size()-1).toString()<<endl;
		if (fragBB(0).distance2(bbAts(0)) > 4){
			continue;
		}
		//MSLOUT.stream() << "RMSD Pre"<<endl;

		double rmsd = fragBB.getAtomPointers().rmsd(bbAts);
		//MSLOUT.stream() << "RMSD: "<<rmsd<<endl;
		if (rmsd > _rmsdTol){
			continue;
		}
		//fragBB.applySavedCoor("pre");

		//MSLOUT.stream() << "RMSD Post"<<endl;
		matchIndex++;
	
		if (_regex != ""){
			if (!boost::regex_search(matchSeq.c_str(), boost::regex(_regex))){
				MSLOUT.stream() << "RegEx NOT Matched. "<<matchSeq<<endl;
				continue;
			} else {
				MSLOUT.stream() << "RegEx Matched. "<<matchSeq<<endl;
			}
			string key = MslTools::stringf("%06d-%s-%1s_%04d%s-%1s_%04d%s",
										   matchIndex,
										   fragDB(i).getSegID().c_str(),
										   fragDB(i).getChainId().c_str(),
										   fragDB(i).getResidueNumber(),
										   fragDB(i).getResidueIcode().c_str(),
										   fragBB(fragBB.size()-1).getChainId().c_str(),
										   fragBB(fragBB.size()-1).getResidueNumber(),
										   fragBB(fragBB.size()-1).getResidueIcode().c_str());
			matchedSequences[key] = matchSeq;
		}

		fprintf(stdout,"(%4s and chain %1s and resi %3d-%3d)  %8.3f \n",
				fragDB[i]->getSegID().c_str(),
				fragDB[i]->getChainId().c_str(),
				fragDB[i]->getResidueNumber(),
				fragBB(fragBB.size()-1).getResidueNumber(),
				rmsd);

		if (pdbDir != ""){
			fragBB.applySavedCoor("pre");
			Atom &at1 = fragBB(0);
			Atom &at2 = fragBB(fragBB.size()-1);
			string allAtomFileName = MslTools::stringf("%s/%s.pdb", pdbDir.c_str(), at1.getSegID().c_str());

			System allAtomSys;
			allAtomSys.readPdb(allAtomFileName);
			for (uint ats = 0; ats < allAtomSys.getAtomPointers().size(); ats++){
				allAtomSys.getAtom(ats).setSegID("");
			}

			if (!tm.rmsdAlignment(fragBB.getAtomPointers(), bbAts, allAtomSys.getAtomPointers())){
				cerr << "ERROR in alignment2: "<<fragBB.size()<<" "<<bbAts.size()<<endl;
				continue;
			}

			lastResults.push_back(new AtomContainer());
			if (includeFullFile){
				lastResults.back()->addAtoms(allAtomSys.getAtomPointers());
			} else {
				AtomSelection sel2(allAtomSys.getAtomPointers());
				ss.str("");
				char tmpstr2[100];
				sprintf(tmpstr2,"chain %1s and resi %d-%-d", at1.getChainId().c_str(), at1.getResidueNumber(), at2.getResidueNumber());
				ss << tmpstr2;
				AtomPointerVector allAts = sel2.select(ss.str());
				lastResults.back()->addAtoms(allAts);
			}

		} else {
			// Add CA only atoms..
			//lastResults.push_back(new AtomContainer(fragBB));

			// For each position in bbAts add the sequence from matchedSequence[i] to a string
			if (bbAts.size() != matchSeq.length()) {
				cerr << "ERROR 2342 PDBFragments::searchForMatchingFragmentsLinear, bbAts.size() != matchedSeq.length()"<<endl;
			} else {
				for (uint x = 0; x < bbAts.size(); x++){
					string key = MslTools::stringf("%s-%d-%s", bbAts[x]->getChainId().c_str(), bbAts[x]->getResidueNumber(), bbAts[x]->getResidueIcode().c_str());
					//cout<< "Key: "<<key<<" "<<matchSeq.substr(x, 1)<<endl;
					//cout << "\tadding to :"<<matchedSequences[MslTools::stringf("%s-%d-%s", bbAts[x]->getChainId().c_str(), bbAts[x]->getResidueNumber(), bbAts[x]->getResidueIcode().c_str())]<<endl;
					matchedSequences[MslTools::stringf("%s-%d-%s", bbAts[x]->getChainId().c_str(), bbAts[x]->getResidueNumber(), bbAts[x]->getResidueIcode().c_str())] += matchSeq.substr(x, 1);
				}
			}    
		}
		
		fragBB.removeAllAtoms();

		numFrags++;
	}

	return numFrags;
}

/*
Input:    
    _sys                    :   MSL System
    _stemResidues           :   id of "stem" residues inside chain   (current needs to be size of 4, 2 on each side of the locally sampling fragment)
    _numResiduesInFragment :*  number of residues wanted between "stems" (default is calculated from firstOfLastTwoStems.getResidueNumber() - lastOfFirstTwoStems.getResidueNumber() )

Output:
    int number of fragments found.
 */
int PDBFragments::searchForMatchingFragmentsStems(System &_sys, vector<string> &_stemResidues,int _numResiduesInFragment, string _regex, double _rmsdTol){

	if (fragDB.size() <= 4){
		cerr << "ERROR 1232 PDBFragments::searchForMatchingFragments(), fragment database is too small, maybe you forgot to load it using loadFragmentDatabase() ?"<<endl;
		exit(1232);
	}
	int numFrags = 0;
	int illegalQuads = 0;

	// Remove last set of results
	lastResults.clear();

	// Clear our set of matched sequences
	matchedSequences.clear();


	// If stem residues is 0, the use all the atoms of _ats to search
	if (_stemResidues.size() == 2){
	} else {

		// Get two stem vectors
		if (_stemResidues.size() % 2 != 0){
			cerr << "ERROR 1967 PDBFragments::searchForMatchingFragments() , number of stem residues is not even\n";
			exit(1967);
		}


		// Create 2 lists (AtomPointerVectors) of the cterm and nterm stem.
		int numStemResidues =  _stemResidues.size() / 2;
		AtomPointerVector stem1;
		AtomPointerVector stem2;
		for (uint i = 0; i < _stemResidues.size();i++){
			if (i < numStemResidues){
				stem1.push_back(&_sys.getPosition(_stemResidues[i]).getAtom("CA"));
			} else {
				stem2.push_back(&_sys.getPosition(_stemResidues[i]).getAtom("CA"));
			}

		}

		AtomPointerVector stemBBats;
		stemBBats.push_back(&_sys.getPosition(_stemResidues[0]).getAtom("N"));
		stemBBats.push_back(&_sys.getPosition(_stemResidues[0]).getAtom("CA"));
		stemBBats.push_back(&_sys.getPosition(_stemResidues[0]).getAtom("C"));
		stemBBats.push_back(&_sys.getPosition(_stemResidues[_stemResidues.size()-1]).getAtom("N"));
		stemBBats.push_back(&_sys.getPosition(_stemResidues[_stemResidues.size()-1]).getAtom("CA"));
		stemBBats.push_back(&_sys.getPosition(_stemResidues[_stemResidues.size()-1]).getAtom("C"));

		//MSLOUT.stream() << "STEM sizes: "<<stem1.size()<<","<<stem2.size()<<endl;

		// Check for small stems
		if (stem1.size() <= 1 || stem2.size() <= 1){
			cerr<< "PDBFragments::searchForMatchingFragments() . Stems are too small."<<stem1.size()<<","<<stem2.size()<<endl;
		}

		// Use natural sequence length if num residues in fragment not set
		if (_numResiduesInFragment == -1){
			_numResiduesInFragment = stem2(0).getResidueNumber() - stem1(stem1.size()-1).getResidueNumber() - 1;
		}

		MSLOUT.stream() << "Number of residues between stems: "<<_numResiduesInFragment<<endl;

		// Store stem-to-stem distance-squared vector (for filtering candidates below)
		AtomPointerVector stems = stem1 + stem2;
		vector<double> stemDistanceSq;
		for (uint c = 0; c < stem1.size();c++){
			  fprintf(stdout,"C: %1s %4d %3s\n",
			  stem1[c]->getChainId().c_str(),
			  stem1[c]->getResidueNumber(),
			  stem1[c]->getResidueName().c_str());
			for (uint n = 0; n < stem2.size();n++){

				double distSq = stem1(c).distance2(stem2(n));			
				  fprintf(stdout,"\tN: %1s %4d %3s = %8.3f\n",
				  stem2[n]->getChainId().c_str(),
				  stem2[n]->getResidueNumber(),
				  stem2[n]->getResidueName().c_str(),distSq);

				stemDistanceSq.push_back(distSq);
			}
		}


		// BBQTable is used for adding backbone-atoms to a c-alpha trace.
		BBQTable bbqT;
		//bbqT.setDebugFlag(true);
		if (fragType == caOnly){

			// BBQ Table for adding backbone atoms
			bbqT.openReader(bbqTable);
		}

		MSLOUT.stream() << "FragDB.size(): "<<fragDB.size()<<endl;
		// Now loop over all fragments in database checking for ones of the correct size
		double tol = 64; // Tolerance of distance to be deviant from stems in Angstroms^2
		int matchIndex = 0; // Index for keeping track of matches
		for (uint i = 0 ; i < fragDB.size()-(_numResiduesInFragment+stem1.size()+stem2.size());i++){

			// Get proposed ctermStem
			AtomPointerVector ctermStem;
			for (uint n = 0; n < stem1.size();n++){
				ctermStem.push_back(fragDB[i+n]);
			}

			// Get proposed ntermStem
			AtomPointerVector ntermStem;
			for (uint n = 0; n < stem2.size();n++){
				ntermStem.push_back(fragDB[i+stem1.size()+_numResiduesInFragment+n]);
			}

			// Check for same PDB
			if (ctermStem(0).getSegID() != ntermStem(ntermStem.size()-1).getSegID()){
			  //fprintf(stdout,"New PDB: %4s to %4s\n",ctermStem[0]->getSegID().c_str(),ntermStem(ntermStem.size()-1).getSegID().c_str());

				// Jump ahead to move past all pdb1 v pdb2 tests
				i = i + _numResiduesInFragment+stem1.size();
				continue;
			}


			// Check for same Chain
			if (ctermStem(0).getChainId() != ntermStem(ntermStem.size()-1).getChainId()){
				/*
				fprintf(stdout,"New Chain: %4s %1s to %4s %1s\n",
					ctermStem[0]->getSegID().c_str(),
					ctermStem[0]->getChainId().c_str(),					
					ntermStem(ntermStem.size()-1).getSegID().c_str(),
					ntermStem(ntermStem.size()-1).getChainId().c_str());
				*/
				// Jump ahead to move past all chain1 v chain2 tests
				i = i + _numResiduesInFragment+stem1.size();
				continue;
			}


			// Check for a gap (PDB Structures have gaps).
			if ( abs(ctermStem[ctermStem.size()-1]->getResidueNumber() - ntermStem[0]->getResidueNumber()) != _numResiduesInFragment+1){
				/*
				fprintf(stdout," Gap found at residue %4s %1s %4d %3s and %4s %1s %4d %3s\n",
					ctermStem[ctermStem.size()-1]->getSegID().c_str(),
					ctermStem[ctermStem.size()-1]->getChainId().c_str(),
					ctermStem[ctermStem.size()-1]->getResidueNumber(),
					ctermStem[ctermStem.size()-1]->getResidueName().c_str(),
					ntermStem[0]->getSegID().c_str(),
					ntermStem[0]->getChainId().c_str(),
					ntermStem[0]->getResidueNumber(),
					ntermStem[0]->getResidueName().c_str());
				*/
				// Jump ahead to move past gap
				i = i + _numResiduesInFragment+stem1.size();
				continue;
			}
			

			// Distance Filter... are the proposed stems close enough ?
			bool passDistanceFilter = true;
			int index = 0;

			// Check against compiled list of distance-squared's
			for (uint c = 0; c < ctermStem.size();c++){

			  /*
				  fprintf(stdout," Checking residue %1s %4d %3s\n",
				  ctermStem[c]->getChainId().c_str(),
				  ctermStem[c]->getResidueNumber(),
				  ctermStem[c]->getResidueName().c_str());
			  */
				for (uint n = 0; n < ntermStem.size();n++){

					double distSq = ctermStem[c]->distance2(*ntermStem[n]);

					/*
					  fprintf(stdout," \t%1s %4d %3s = %8.3f vs %8.3f\n",
					  ntermStem[n]->getChainId().c_str(),
					  ntermStem[n]->getResidueNumber(),
					  ntermStem[n]->getResidueName().c_str(),distSq,stemDistanceSq[index]);
					*/
					if (fabs(stemDistanceSq[index++] - distSq) > tol){
						passDistanceFilter = false;
						break;
					}
				}

				if (!passDistanceFilter){
					break;
				}
			}



			// Continue if the distance filter not passed
			if (!passDistanceFilter) {
				continue;
			}
			MSLOUT.stream() << "PASSED DISTANCE FILTER!"<<endl;

			// align and print winning fragment

			// Make a copy of frag stem, so not to effect original atoms
			AtomPointerVector fragStem;
			for (uint ct = 0; ct < ctermStem.size();ct++){
				fragStem.push_back(new Atom(ctermStem(ct)));
			}
			for (uint nt = 0; nt < ntermStem.size();nt++){
				fragStem.push_back(new Atom(ntermStem(nt)));
			}

			// Align the fragStem to stem to see if it matches well enough...
			Transforms tm;
			fragStem.saveCoor("pre");
			tm.rmsdAlignment(fragStem,stems);
			
			double rmsd = fragStem.rmsd(stems);
			MSLOUT.stream() << "RMSD: "<<rmsd<<endl;

			// Continue if the RMSD filter not passed
			if (rmsd > _rmsdTol){
				continue;
			}



			// Adding clash filter?
			//if (_checkForClashes){
			  
			//}



			matchIndex++;
			fragStem.applySavedCoor("pre");

			//lastResults->writePdb("/tmp/preAdd.pdb");

			// Create new atoms for middle atoms..
			AtomPointerVector tmp;
			index = 1;
			string matchSeq = "";
			for (uint f = 0; f < ctermStem.size();f++){
			        //fragStem(f).setResidueNumber(index++);
				//fragStem(f).setResidueName("FRG");
				//fragStem(f).setChainId("A");
				//fragStem(f).setSegID("");
				tmp.push_back(&fragStem(f));
				//cout << "C-ADD: "<<fragStem(f)<<" "<<fragStem(f).getSegID()<<endl;
				matchSeq += MslTools::getOneLetterCode(fragStem(f).getResidueName());
			}


			for (int f = 0; f < _numResiduesInFragment;f++){
			        Atom *a =new Atom(fragDB(i+ctermStem.size()+f));
				//a->setCoor(fragDB(i+ctermStem.size()+f).getCoor());
				//a->setResidueNumber(index++);
				//a->setResidueName("FRG");
				//a->setChainId("A");
				//a->setSegID("");
				tmp.push_back(a);
				//cout << "F-ADD: "<<*a<<" "<<a->getSegID()<<endl;

				matchSeq += MslTools::getOneLetterCode(fragDB(i+ctermStem.size()+f).getResidueName());
			} 		


			for (uint f = 0; f < ntermStem.size();f++){
				//fragStem(f+ctermStem.size()).setResidueNumber(index++);
				//fragStem(f+ctermStem.size()).setResidueName("FRG");
				//fragStem(f+ctermStem.size()).setChainId("A");
				//fragStem(f+ctermStem.size()).setSegID("");
				tmp.push_back(&fragStem(f+ctermStem.size()));
				matchSeq += MslTools::getOneLetterCode(fragStem(f+ctermStem.size()).getResidueName());
				//cout << "N-ADD: "<<fragStem(f+ctermStem.size())<<" "<<fragStem(f+ctermStem.size()).getSegID()<<endl;
			}
			string key = MslTools::stringf("%06d-%s-%1s_%04d%s-%1s_%04d%s",
						       matchIndex,
						       fragDB(i+ctermStem.size()).getSegID().c_str(),
						       fragDB(i+ctermStem.size()).getChainId().c_str(),
						       fragDB(i+ctermStem.size()).getResidueNumber(),
						       fragDB(i+ctermStem.size()).getResidueIcode().c_str(),
						       fragDB(i+ctermStem.size()+_numResiduesInFragment-1).getChainId().c_str(),
						       fragDB(i+ctermStem.size()+_numResiduesInFragment-1).getResidueNumber(),
						       fragDB(i+ctermStem.size()+_numResiduesInFragment-1).getResidueIcode().c_str());

			matchedSequences[key] = matchSeq;

			if (_regex != ""){


			  if (!boost::regex_search(matchSeq.c_str(),boost::regex(_regex))){
			    MSLOUT.stream() << "RegEx NOT Matched. "<<matchSeq<<endl;
			    continue;
			  } else {
			    MSLOUT.stream() << "RegEx Matched."<<endl;
			  }

			}
			

			Chain tmpChain;
			tmpChain.addAtoms(tmp);			

			
			tm.rmsdAlignment(fragStem,stems,tmpChain.getAtomPointers());


			fprintf(stdout,"(%4s and chain %1s and resi %3d-%3d)  %8.3f ",
				ctermStem[0]->getSegID().c_str(),
				ctermStem[0]->getChainId().c_str(),
				ctermStem[0]->getResidueNumber(),
				ntermStem[ntermStem.size()-1]->getResidueNumber(),
			        rmsd);

			
			bool successful = true;
			if (fragType == caOnly){
			  if (pdbDir != ""){
			    
			              Atom &at1 = tmpChain.getAtom(0);
			              Atom &at2 = tmpChain.getAtom(tmpChain.atomSize()-1);

				      string allAtomFileName = MslTools::stringf("%s/%s.pdb",pdbDir.c_str(),at1.getSegID().c_str());
				      MSLOUT.stream() << "Opening "<<allAtomFileName<<endl;
				      System allAtomSys;
				      allAtomSys.readPdb(allAtomFileName);
				      for (uint ats = 0; ats < allAtomSys.getAtomPointers().size();ats++){
					allAtomSys.getAtom(ats).setSegID("");
				      }

				      AtomSelection sel2(allAtomSys.getAtomPointers());

			    	      stringstream ss;
				      char tmpstr[100];
				      //sprintf(tmpstr,"chain %1s and resi %d-%-d and name CA",at1.getChainId().c_str(),at1.getResidueNumber(),at2.getResidueNumber());
				      sprintf(tmpstr,"chain %1s and resi %d-%-d and name CA",	ctermStem[0]->getChainId().c_str(),ctermStem[0]->getResidueNumber(),ntermStem[ntermStem.size()-1]->getResidueNumber());
				      ss << tmpstr;
				      AtomPointerVector caAts = sel2.select(ss.str());


				      if (!tm.rmsdAlignment(caAts,tmpChain.getAtomPointers(),allAtomSys.getAtomPointers())){
					MSLOUT.stream() << "Problem aligning all atoms using the C-alpha trace"<<endl;
					MSLOUT.stream() << "PDB: "<<allAtomFileName<<endl;
					MSLOUT.stream() << "\tTrying to align with: "<<ss.str()<<endl;
					MSLOUT.stream() << "\tSystem: "<<allAtomSys.getSizes();
					MSLOUT.stream() << "\tSelected: "<<caAts.size()<<" atoms using '"<<ss.str()<<"'"<<endl;
					MSLOUT.stream() << "\tReference: "<<tmpChain.getAtomPointers().size()<<" atoms"<<endl;
					MSLOUT.stream() << tmpChain.getAtomPointers();
					successful = false;
					continue;
				      } 

				      // Get BB atoms from stem-equivalent residues
				      AtomPointerVector allAts_stemBBats;
				      //Position &stem1eq  = allAtomSys.getPosition(MslTools::getPositionId(at1.getChainId(),at1.getResidueNumber(),at1.getResidueIcode()));
				      Position &stem1eq  = allAtomSys.getPosition(MslTools::getPositionId(ctermStem[0]->getChainId(),ctermStem[0]->getResidueNumber(),ctermStem[0]->getResidueIcode()));
				      allAts_stemBBats.push_back(&stem1eq.getAtom("N"));
				      allAts_stemBBats.push_back(&stem1eq.getAtom("CA"));
				      allAts_stemBBats.push_back(&stem1eq.getAtom("C"));

				      //Position &stem2eq  = allAtomSys.getPosition(MslTools::getPositionId(at2.getChainId(),at2.getResidueNumber()+stem2.size()-1,at2.getResidueIcode()));
				      Position &stem2eq  = allAtomSys.getPosition(MslTools::getPositionId(ntermStem[ntermStem.size()-1]->getChainId(),ntermStem[ntermStem.size()-1]->getResidueNumber(),ntermStem[ntermStem.size()-1]->getResidueIcode()));
				      allAts_stemBBats.push_back(&stem2eq.getAtom("N"));
				      allAts_stemBBats.push_back(&stem2eq.getAtom("CA"));
				      allAts_stemBBats.push_back(&stem2eq.getAtom("C"));
				      	
				      double bbRMSD = allAts_stemBBats.rmsd(stemBBats);

				      MSLOUT.stream() << "BB RMSD: "<<bbRMSD<<endl;
				      fprintf(stdout, "%8.3f",bbRMSD);
				      //if (bbRMSD > 1.51){
				      //if (bbRMSD > 0.3){
				      if (bbRMSD > _rmsdTol*3){
					successful = false;
					continue;
				      }
				      MSLOUT.stream() << "ADDDDDDDDDDDDDDDDDING COORDS"<<endl;


				      lastResults.push_back(new AtomContainer());
				      if (includeFullFile){
					lastResults.back()->addAtoms(allAtomSys.getAtomPointers());
				      } else {
					ss.str("");
					char tmpstr2[100];
					sprintf(tmpstr2,"chain %1s and resi %d-%-d",at1.getChainId().c_str(),at1.getResidueNumber(),at2.getResidueNumber());
					ss << tmpstr2;
					AtomPointerVector allAts = sel2.select(ss.str());
					lastResults.back()->addAtoms(allAts);
				      }






			  } else if (bbqTable != ""){
			  
				illegalQuads += bbqT.fillInMissingBBAtoms(tmpChain);

				if (illegalQuads > 0){
				  successful = false;
			  	  exit(0);
				}

				lastResults.push_back(new AtomContainer(tmpChain.getAtomPointers()));
			  }
				
			} else {
			  MSLOUT.stream() <<" ALL ATOMS "<<endl;
			}
			fprintf(stdout,"\n");

			fragStem.applySavedCoor("pre");

			if (successful){
				numFrags++;
			}

			/*					
			stringstream ss;
			PDBWriter pout;
			pout.open(ss);
			pout.write(sys.getAtomPointers());
			pout.close();

			matchingFragments.push(pair<double,string>(rmsd,ss.str()));
			*/


		}

		//fprintf(stdout, "Best RMSD: %8.3f\n",matchingFragments.top().first);


		fprintf(stdout,"Number of succesful fragments: %10d, illegal quads: %10d\n",numFrags,illegalQuads);

	}
	
	
	return numFrags;
}


