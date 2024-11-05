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
// MSL Includes
#include "System.h"
#include "MslTools.h"
#include "AtomSelection.h"
#include "AtomBondBuilder.h"
#include "Transforms.h"
#include "OptionParser.h"
#include "MslOut.h"
#include "ChiStatistics.h"
#include "RandomNumberGenerator.h"
#include "Atom3DGrid.h"
#include "glycanTransfer.h"


// STL Includes
#include<iostream>
#include<vector>
#include<map>

using namespace std;

using namespace MSL;


// MslOut 
static MslOut MSLOUT("glycanTransfer");



// In options: 	vector<DoF> edits; // atom selection one

int main(int argc, char *argv[]) {	


	// Option Parser
	Options opt = setupOptions(argc,argv);
	
	System *pdb = new System();
	pdb->readPdb(opt.pdb);

	AtomContainer pdbCon(pdb->getAtomPointers());
	
	// Create a Atom3DGrid
	AtomSelection bb_sel(pdbCon.getAtomPointers());
	
	Atom3DGrid grid(bb_sel.select("name CA+CB+N+C+O"),3.0);

	// ChiStatistics
	ChiStatistics chi;
	chi.read(opt.doffile);

	// Extract new residue from newRes pdb file
	System newResidue;
	newResidue.readPdb(opt.newRes);
	newResidue.getAtomPointers().saveCoor("init");

	// Build bonds
	AtomBondBuilder abb;
	abb.buildConnections(newResidue.getAtomPointers());	

	// Parse alignable and non-alignable atoms
	AtomSelection sel(newResidue.getAtomPointers());
	AtomPointerVector newResAlignAts = sel.select(opt.alignAtoms);
	AtomPointerVector newResNotAlignAts = sel.select("not "+opt.alignAtoms);
	AtomPointerVector newResRMSDAts;
	AtomContainer newResRMSDCheck;
	if (opt.rmsdAtoms != ""){
	  newResRMSDAts = sel.select(opt.rmsdAtoms);
	  newResRMSDCheck.addAtoms(newResRMSDAts);
	  cout << "NEW RES RMSD CHECK: "<<endl<<newResRMSDCheck.toString()<<endl;
	}

	for (uint i =0; i < newResidue.atomSize();i++){
	  cout << "ATom name: "<<newResidue.getAtom(i).getName()<<"."<<endl;
	}

	// For each position
	for (uint i = 0; i < opt.positions.size();i++){

	  if (!pdb->positionExists(opt.positions[i])){
	    cerr << "ERROR 2222 position "<<opt.positions[i]<< " does not exist in "<<opt.pdb<<" skipping this position"<<endl;
	    continue;
	  }
	  cout << "Working on position: "<<opt.positions[i]<<endl;

	  grid.resetAllHiddenFlags();
	  grid.setPositionToSkip(pdb->getPosition(opt.positions[i]).getPositionId());
	  
	  // Select atoms to align
	  AtomSelection sel(pdb->getPosition(opt.positions[i]).getAtomPointers());
	  AtomPointerVector posAlignAts = sel.select(opt.alignAtoms);
	  if (newResAlignAts.size() != posAlignAts.size()){
	    cerr << "ERROR 333 position "<<opt.positions[i]<< " has "<<posAlignAts.size()<< " atoms and residue to add has "<<newResAlignAts.size()<<" skipping this position" <<endl;
	    cout << "Position Align Ats" <<endl;
	    cout << posAlignAts.toString()<<endl;
	    cout << "New Res Align Ats" <<endl;
	    cout << newResAlignAts.toString()<<endl;
	    continue;
	  }

	  /* 

	     Search for best conformation for the newResidue on the pdb at position[i]
	     
	      opt.maxIter   = X ;  [ X = number of iterations to sample each glycan ]
	     
              opt.randomize = 0 ;  [ sample dihedrals in order from Asn ]
	      opt.randomize = 1 ;  [ sample dihedrals in random order ] * looks like first iteration that it does not do random...
	      opt.randomize = 2 ;  [ opt.randomize = 1 AND up-to the last 1/4 of iterations and past the first 2, re-init structure , else use globalMin ] 
	      opt.randomize = 3 ;  [ start by randomizing all angles one at a time, then go from top to bottom ]
	     
	   */
	  
	  int globalMinClashes = MslTools::intMax;
	  double globalMinRMSD = MslTools::doubleMax;

	  RandomNumberGenerator rng;
	  cout << "RNG SEED: "<<rng.getSeed()<<endl;

	  for (uint iter = 0; iter < opt.maxIter; iter++){
	    if (globalMinClashes == 0) continue;

	    
	    double minRMSD = MslTools::doubleMax;
	    int minClashes = MslTools::intMax;
	    if (opt.randomize == 2){

	      // Up-to the last 1/4 of iterations and past the first 2 iterations , re-init structure to try again.
	      if (iter < opt.maxIter - int(opt.maxIter/4)){
		if (iter > 1){
		  newResidue.getAtomPointers().applySavedCoor("init");
		}
	      } else {
		newResidue.getAtomPointers().applySavedCoor("globalMin");
		minRMSD = globalMinRMSD;
		minClashes = globalMinClashes;
	      }

	    }


	    // For each flexible dihedral..
	    map<int,bool> dihedral_used;
	    for (uint d = 0; d < opt.DoFs.size();d++){
	      if (globalMinClashes == 0) continue;
	      
	      int d_index = d;


	      // If randomize pick one.
	      if (opt.randomize >= 1 && iter > 0){
		map<int,bool>::iterator it;
		int tries = 0;
		int new_d = 0;
		do {
		  new_d = rng.getRandomInt(opt.DoFs.size()-1);
		  it = dihedral_used.find(new_d);
		  tries++;
		}while (tries < 100000 && it != dihedral_used.end());

		if (tries == 100000) break;
		dihedral_used[new_d] = true;

		d_index = new_d;
	      }


	      fprintf(stdout, ">>>>>Sampling position %10s dihedral %10d %10d %8.2f %8.2f %8.2f %10d %10d\n",opt.positions[i].c_str(),d,d_index,opt.DoFs[d_index].range_start,opt.DoFs[d_index].range_stop,opt.DoFs[d_index].range_step, minClashes, globalMinClashes);

	      // Range of dihedrals
	      for (double r = opt.DoFs[d_index].range_start; r < opt.DoFs[d_index].range_stop; r+=opt.DoFs[d_index].range_step){
		if (globalMinClashes == 0) continue;

		
		// Align to pdb position[i]
		Transforms t;
		if (!t.rmsdAlignment(newResAlignAts, posAlignAts, newResidue.getAtomPointers())){
		  cerr << "ERROR 4444 alignment error with position "<<opt.positions[i]<<" skipping transformation "<<endl;
		  continue;
		}
		newResidue.writePdb("newResidue1.pdb");
		// Create conformation
		setDihedral(newResidue,opt.DoFs[10],r);
		cout << "R: "<<r<<endl;
		newResidue.writePdb("newResidue2.pdb");
		exit(0);
		// Clash check
		int numClashes = countClashes(newResNotAlignAts, grid, false, opt.clashTol*100);
		//int numClashes = countClashes(newResNotAlignAts,pdb->getAtomPointers(),false);
		//newResidue.writePdb("newResidue.pdb");
		

		if (numClashes < opt.clashTol) {
		  
		  if (opt.rmsdAtoms != ""){
		    double rmsd = newResRMSDCheck.getAtomPointers().rmsd(newResRMSDAts);

		    if (rmsd < minRMSD){
		      //cout << "MinRMSD = "<<rmsd<<" clashes: "<<numClashes<<endl;

		      fprintf(stdout, "MinRMSD found: %8.2f for dihedral %03d using angle %8.2f chis [ ",rmsd,d,r);
		      for (uint c = 1; c <= chi.getNumberChis(newResidue.getResidue(0)); c++){
			double chi_angle = chi.getChi(newResidue.getResidue(0), c);
			if (chi_angle == MslTools::doubleMax){
			  chi_angle = 0.0;
			}
			fprintf(stdout," %6.1f",chi_angle);
		      }
		      fprintf(stdout,"\n");

		      newResidue.getAtomPointers().saveCoor("min");
		      minRMSD = rmsd;

		      if (minRMSD < globalMinRMSD){
			newResidue.getAtomPointers().saveCoor("globalMin");
			globalMinRMSD = minRMSD;
		      }

		    } else {
		      //cout << "RMSD is "<<rmsd<<endl;
		    }
		  } // END IF opt.rmsdAtoms 
		  else {
		  } // ELSE (opt.rmsdAtoms == "")


		} else { // ELSE numClashes >= opt.clashTol
		  //MSLOUT.stream() << "TOO MANY BUMPS: "<<numClashes<<endl;
		}

		if (numClashes < minClashes){
		      cout << "MinClashes: "<<numClashes<<endl;
		      newResidue.getAtomPointers().saveCoor("min");
		      minClashes = numClashes;
		} // IF new minClashes
		
		// Best model is minimum clashes...
		if (minClashes < globalMinClashes){
		      cout << "GlobalMinClashes: "<<minClashes<<endl;
		      newResidue.getAtomPointers().saveCoor("globalMin");
		      globalMinClashes = minClashes;
		}		

		
	      } // END FOR DOF range loop
	      fprintf(stdout, "\tFull conformation sampled.\n");
	      newResidue.getAtomPointers().applySavedCoor("min");

	    } // FOR DIHEDRALS  
	  } // FOR ITER
	  
	  // Retreive best solution
	  newResidue.getAtomPointers().applySavedCoor("globalMin");
	  cout << "Global minimum clashes: "<<globalMinClashes<<endl;

	  // Create a new pdb with this conformation
	  int index = pdb->getPositionIndex(opt.positions[i]);
	  System *newPdb = new System();
	  for (uint p = 0; p < index;p++){
	    newPdb->addAtoms(pdb->getPosition(p).getAtomPointers());
	  }
	  // get chain and residue number from position i
	  string chain = pdb->getPosition(opt.positions[i]).getChainId();
	  int resn     = pdb->getPosition(opt.positions[i]).getResidueNumber();
	  string icode = pdb->getPosition(opt.positions[i]).getResidueIcode();
	  AtomPointerVector newAts;
	  for (uint a = 0; a < newResidue.getAtomPointers().size();a++){
	    //string id = MslTools::stringf("%s,%d,%s,%s", chain.c_str(),resn,newResidue.getAtom(a).getResidueName().c_str(),newResidue.getAtom(a).getName().c_str());
	    string id = MslTools::getAtomId(chain,resn,icode,newResidue.getAtom(a).getName());
	    Atom *at = new Atom(id);
	    at->setCoor(newResidue.getAtom(a).getCoor());
	    at->setResidueName(newResidue.getAtom(a).getResidueName());
	    //cout << "ADDING "<<id<<endl;
	    newPdb->addAtom(*at);
	    newAts.push_back(at);
	  }
	  
	  
	  for (uint p = index+1;p < pdb->positionSize();p++){
	    newPdb->addAtoms(pdb->getPosition(p).getAtomPointers());
	  }
	  newPdb->writePdb(MslTools::stringf("%s_position_%04d.pdb",MslTools::getFileName(opt.outPdb).c_str(), pdb->getPosition(opt.positions[i]).getResidueNumber()));
	  		   
	  // Add new atoms to grid 
	  grid.addAtoms(newAts);
	  
	  // Make 'pdb' this newPdb...
	  delete(pdb);
	  pdb = newPdb;
	  newPdb = NULL;
	  newResidue.getAtomPointers().clearSavedCoor("min");
	  newResidue.getAtomPointers().clearSavedCoor("globalMin");

	  
	} // FOR Positions(i)

	pdb->writePdb(opt.outPdb);

	
}

Options setupOptions(int theArgc, char * theArgv[]){
	Options opt;

	OptionParser OP;


	OP.setRequired(opt.required);
	OP.setAllowed(opt.optional);
	OP.setDefaultArguments(opt.defaultArgs);
	OP.autoExtendOptions();
	
	OP.readArgv(theArgc, theArgv);
	opt.configfile = OP.getString("configfile");
	if (opt.configfile != "") {
		OP.readFile(opt.configfile);
		if (OP.fail()) {
		  cerr <<  "ERRROR 1111 Cannot read configuration file " << opt.configfile <<endl;
			exit(1);
		}
	}

	if (OP.countOptions() == 0){
		cout << "Usage:" << endl;
		cout << endl;
		exit(0);
	}

	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}


	opt.newRes = OP.getString("newRes");
	if (OP.fail()){
		cerr << "ERROR 1111 newRes not specified.\n";
		exit(1111);
	}

	opt.positions = OP.getMultiString("positions");
	if (OP.fail()){
		cerr << "ERROR 1111 positions not specified.\n";
		exit(1111);
	}

	opt.doffile = OP.getString("doffile");
	DegreeOfFreedomReader chi;
	if  (OP.fail()){
	  opt.doffile = "";
	} else {
	  chi.read(opt.doffile);
	}
	vector<string> ranges = OP.getMultiString("range");

	// Convert dihedrals and ranges into DoF structs.
	int index = 0;
	while (true) {
		DoF tmp;
		tmp.type = "dihedral";

		tmp.atomNames = OP.getStringVector("dihedrals", index);
		if (OP.fail()) {
			break;
		}

		// Parse out atom names from pre-defined chi angles  RESNAME,CHIN (e.g. ARG,chi4)
		if (opt.doffile != ""  && tmp.atomNames.size() == 1){
		  
		  vector<string> items = MslTools::tokenize(tmp.atomNames[0],",");
		  tmp.atomNames = chi.getSingleDegreeOfFreedom(items[0],items[1]);
		  if (tmp.atomNames.size() == 0){
		    cerr << "ERROR 1111 could not find DoF '"<<items[0]<<","<<items[1]<<"' in "<<opt.doffile<<endl;
		    exit(1111);
		  }

		}

		if (tmp.atomNames.size() != 4) {
		  cerr << "ERROR 1111 Didn't get 4 atom names for a dihedral specification, pre-defined names require '--doffile FILE'"<<endl;
		  exit(1111);
		}

		if (ranges.size() != 0){
		  if (!OP.fail()){
		    vector<string> items = MslTools::tokenize(ranges[index], ":");
		    if (items.size() != 3){
		      cerr << "ERROR range for dihedrals is not valid, does not have 3 tokens, delimited by ':'"<<endl;
		      exit(1111);
		    }
		    tmp.range_start = MslTools::toDouble(items[0], "range_start to double");
		    tmp.range_stop  = MslTools::toDouble(items[1], "range_stop to double");
		    tmp.range_step   = MslTools::toDouble(items[2], "range_step to double");
		  }
		}

		opt.DoFs.push_back(tmp);
		index++;
	}


	
	opt.alignAtoms = OP.getString("alignAtoms");
	if (OP.fail()){
	  opt.alignAtoms = "name N+CA+C";
	}

	opt.rmsdAtoms = OP.getString("rmsdAtoms");
	if (OP.fail()){
	  opt.rmsdAtoms = "";
	}
	opt.clashTol = OP.getInt("clashTol");
	if (OP.fail()){
	  opt.clashTol = 10;
	}
	opt.outPdb = OP.getString("outPdb");
	if (OP.fail()){
	  opt.outPdb = "out.pdb";
	}

	opt.maxIter = OP.getInt("maxIter");
	if (OP.fail()){
	  opt.maxIter =1;
	}
	opt.randomize = OP.getInt("randomize");
	if (OP.fail()){
	  opt.randomize = 0;
	}
	return opt;
}



void setDihedral(System &_sys, DoF &_dih, double _value){
  Transforms T;

  if (_dih.atomNames.size() != 4) {
    cerr << "ERROR: Type dihedral/improper requires 4 atom name " <<_dih.atomNames.size()<< endl;
    exit(1);
  }

  
  Atom * pAtom1 = NULL;
  if (_sys.atomExists(_dih.atomNames[0])) {
    pAtom1 = &(_sys.getLastFoundAtom());
  } else {

    // If dihedral is only defined as an atom name and system has only 1 residue allow this to work..
    if (_sys.positionSize() == 1){
      if (_sys.getPosition(0).atomExists(_dih.atomNames[0])){
	pAtom1 = &(_sys.getPosition(0).getLastFoundAtom());
      }
    } else {
      cout << "ERROR: Atom1 " << _dih.atomNames[0] << " not found" << endl;
      exit(1);
    }
  }
  Atom * pAtom2 = NULL;
  if (_sys.atomExists(_dih.atomNames[1])) {
    pAtom2 = &(_sys.getLastFoundAtom());
  } else {

    // If dihedral is only defined as an atom name and system has only 1 residue allow this to work..
    if (_sys.positionSize() == 1){
      if (_sys.getPosition(0).atomExists(_dih.atomNames[1])){
	pAtom2 = &(_sys.getPosition(0).getLastFoundAtom());
      }
    } else {
      cout << "ERROR: Atom2 " << _dih.atomNames[1] << " not found" << endl;
      exit(1);
    }
  }
  Atom * pAtom3 = NULL;
  if (_sys.atomExists(_dih.atomNames[2])) {
    pAtom3 = &(_sys.getLastFoundAtom());
  } else {

    // If dihedral is only defined as an atom name and system has only 1 residue allow this to work..
    if (_sys.positionSize() == 1){
      if (_sys.getPosition(0).atomExists(_dih.atomNames[2])){
	pAtom3 = &(_sys.getPosition(0).getLastFoundAtom());
      }
    } else {

      cout << "ERROR: Atom3 " << _dih.atomNames[2] << " not found" << endl;
      exit(1);
    }
  }
  Atom * pAtom4 = NULL;
  if (_sys.atomExists(_dih.atomNames[3])) {
    pAtom4 = &(_sys.getLastFoundAtom());
  } else {

    // If dihedral is only defined as an atom name and system has only 1 residue allow this to work..
    if (_sys.positionSize() == 1){
      if (_sys.getPosition(0).atomExists(_dih.atomNames[3])){
	pAtom4 = &(_sys.getPosition(0).getLastFoundAtom());
      }
    } else {

      cout << "ERROR: Atom4 " << _dih.atomNames[3] << " not found" << endl;
      exit(1);
    }

  }
  if (_dih.type == "dihedral") {
    T.setDihedral(*pAtom1, *pAtom2, *pAtom3, *pAtom4, _value);
  }

}
int countClashes(AtomPointerVector &_ats, Atom3DGrid &_grid, bool _checkInternal, int _clashTol){
  AtomPointerVector neighbors_on_grid;
  _grid.getNeighbors(_ats,neighbors_on_grid);
  MSLOUT.stream() << "NUM NEIGHBORS: "<<_ats.size()<<" "<<neighbors_on_grid.size()<<endl;
  return countClashes(_ats, neighbors_on_grid,_checkInternal, _clashTol);
}

int countClashes(AtomPointerVector &_ats, AtomPointerVector &_b, bool _checkInternal, int _clashTol){

  int clashes = 0;
  for (uint at2 = 0; at2 < _ats.size();at2++){
  	  
    for (uint at1 = 0; at1 < _b.size();at1++){
      
      if (_b(at1).getResidueName() != "ASX" &&
	  _b(at1).getName() != "CA"  &&
	  _b(at1).getName() != "CB"  &&
	  _b(at1).getName() != "C"  &&
	  _b(at1).getName() != "O"  &&
	  _b(at1).getName() != "N" ) continue;


      if (_b(at1).getSelectionFlag("hidden")) continue;
      
      if (_b(at1).distance2(_ats(at2)) < 9){
	//cout << "Atom "<<_ats(at2) << " and " << _b(at1) << " are this close: "<<_b(at1).distance(_ats(at2))<<endl;
	clashes++;
	if (clashes > _clashTol) {
	  return clashes;
	}
				
      }
    } // FOR AT 1

    if (_checkInternal) {

      // Really close internal clashes as well..
      for (uint at3 = at2+2; at3 < _ats.size();at3++){
	if (_ats(at2).getName().substr(_ats(at2).getName().size()-1,1) == _ats(at3).getName().substr(_ats(at3).getName().size()-1,1)){ 
	  //cout << "Skipping "<<_ats(at2).getName()<<" and "<<_ats(at3).getName()<<endl;
	  continue;
	}
	// Ignore Carbon-Carbon bonds.
	if (_ats(at2).getName().substr(0,1) == "C" || _ats(at3).getName().substr(0,1) == "C"){ 
	  continue;
	}
	if (_ats(at2).distance(_ats(at3)) < 2.5){
	  //cout << "CLASH "<<_ats(at2).getName()<<" and "<<_ats(at3).getName()<<endl;
	  clashes++;
	}	
      } // FOR AT3
    } // IF _checkInternal

  } // FOR AT2

  return clashes;
}
