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
#include "SysEnv.h"
#include "MslTools.h"
#include "OptionParser.h"
#include "MslOut.h"
#include "AtomSelection.h"
#include "AtomBondBuilder.h"
#include "DegreeOfFreedomReader.h"
#include "getDihedrals.h"

// STL Includes
#include<iostream>
#include<map>
#include<string>
#include<vector>


using namespace std;
using namespace MSL;

// MslOut 
static MslOut MSLOUT("getGlycanDihedrals");

int main(int argc, char *argv[]){
	
	// Option Parser
	Options opt = setupOptions(argc,argv);

	// Add pdb to pdblist
	if (opt.pdb != ""){
		opt.pdblist.push_back(opt.pdb);
	}

	
	DegreeOfFreedomReader dof;
	dof.read(opt.dofFile);
	std::map<std::string, std::map<std::string, std::vector<std::string> > > dofvalues = dof.getDegreesOfFreedom();

	for (uint s = 0; s < opt.pdblist.size();s++){
		cout << "Working on " << opt.pdblist[s]<<endl;

		// Read PDB
		System sys;
		sys.readPdb(opt.pdblist[s]);
		
		AtomSelection sel(sys.getAtomPointers());
		AtomPointerVector ats = sel.select(opt.selection);
		System select_sys;
		select_sys.addAtoms(ats);

		AtomBondBuilder abb;
		abb.buildConnections(select_sys.getAtomPointers());


		// Fix glycan atom names here
		//f (!fixAtomNames(select_sys)){
		// cerr << "ERROR 2343 fixing glycan atom names\n";
		// exit(2343);
		//


		string filename = MslTools::getFileName(opt.pdblist[s]);
		map<string,bool> angleMap;
		for (uint i = 0 ; i < select_sys.positionSize();i++){		

		  cout << "For position: "<<i<<endl;
		  Residue & n   = select_sys.getResidue(i);
		  MSLOUT.stream() << "Working on n residue: "<<n.toString()<<endl;



		  // Find directly connected residues.
		  std::vector<int> neighboring_residues = select_sys.getConnectedPositions(*n.getParentPosition());
		  MSLOUT.stream()<<"  Number of neighboring_residues: "<<neighboring_residues.size()<<endl;
		  for (uint j = 0; j < neighboring_residues.size();j++){
		    Residue &m = select_sys.getResidue(neighboring_residues[j]);

		    // Enforce order of residue numbers .. this may cause problems for some files.
		    if (n.getResidueNumber() > m.getResidueNumber()){
		      MSLOUT.stream() << "Skipping "<<n.getResidueNumber()<<" and "<<m.getResidueNumber()<<" because they are out of numerical order"<<endl;
		      continue;
		    }

		    MSLOUT.stream() << "\tWorking on m residue: "<<m.toString()<<endl;
		    string residuePair = MslTools::stringf("%3s-%3s",n.getResidueName().c_str(), m.getResidueName().c_str());
		    std::map<std::string, std::vector<std::string> > angles = dofvalues[residuePair];
		    if (angles.size() == 0){
		      MSLOUT.stream() << "Residue pair: "<<residuePair<<" for positions "<<n.toString()<<" AND "<<m.toString()<< " does not have a DOF definition."<<endl;
		      continue;
		    }

		    std::map<std::string, std::vector<std::string> >::iterator it;
		    for (it = angles.begin();it != angles.end();it++){
		      string angleName = it->first;
		      vector<string> atomNames = it->second;
		      MSLOUT.stream() << "Working on angle: "<<angleName<<endl;		      



		      AtomPointerVector atomPtrs;
		      for (uint at = 0; at < atomNames.size();at++){

			// See which residue this atom is part of
			if (atomNames[at][0] != '+'){
			  if (!n.atomExists(atomNames[at])){
			    MSLOUT.stream() <<"ERROR 3333 AtomName : "<<atomNames[at]<<" does not exist in residue: "<<n.toString()<<endl;
			    continue;
			  } else {
			      atomPtrs.push_back(&n.getLastFoundAtom());
			    }
			} else {
			  if (!m.atomExists(atomNames[at].substr(1))){
			    MSLOUT.stream() <<"ERROR 3334 AtomName : "<<atomNames[at]<<" does not exist in residue: "<<m.toString()<<endl;
			    continue;
			  } else {
			      atomPtrs.push_back(&m.getLastFoundAtom());
			  }

			}

		      }

		      if (atomPtrs.size() != 4){ 
			MSLOUT.stream() << "Could not get 4 atom pointers from atomNames: "<<n.toString()<<" and "<<m.toString()<<endl;
			continue;
			//exit(23423);
		      }
		      
		      // Check that these atoms are covalently bonded.
		      if (!(atomPtrs[0]->isBoundTo(atomPtrs[1]) &&
			    atomPtrs[1]->isBoundTo(atomPtrs[2]) &&
			    atomPtrs[2]->isBoundTo(atomPtrs[3])
			    )
			 ){
			MSLOUT.stream() << "Angle: "<<angleName<<" not applicable, some of these atoms are not bonded: "<<atomPtrs[0]->getName()<<" "<<atomPtrs[1]->getName()<<" "<<atomPtrs[2]->getName()<<" "<<atomPtrs[3]->getName()<<endl;
			continue;
		      }

		      // Finally get the dihedral
		      double dih = atomPtrs[0]->dihedral(*atomPtrs[1], *atomPtrs[2], *atomPtrs[3]);




		      fprintf(stdout, "%s %1s %3s %3d%1s %3s %3d%1s %10s %6.2f\n",filename.c_str(),n.getChainId().c_str(),
			      n.getResidueName().c_str(),n.getResidueNumber(),n.getResidueIcode().c_str(),
			      m.getResidueName().c_str(),m.getResidueNumber(),m.getResidueIcode().c_str(),
			      angleName.c_str(),
			      dih
			      );		      
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
	OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"
	OP.readArgv(theArgc, theArgv);

	if (OP.countOptions() == 0){
		cout << "Usage:" << endl;
		cout << endl;
		cout << "getGlycanDihedrals --pdb PDB --doffile <DEGREE_OF_FREEDOM_FILE>\n";
		exit(0);
	}

	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		opt.pdblist = OP.getStringVector("pdblist");
		if (OP.fail()){
			cerr << "ERROR 1111 pdb or pdblist not specified.\n";
			exit(1111);
		}
	} else {
		opt.pdblist = OP.getStringVector("pdblist");
	}

	opt.dofFile = OP.getString("doffile");
	if (OP.fail()){
	        SysEnv env;
		opt.dofFile = env.getEnv("MSL_PDB_2_3_DOF");
		cerr << "WARNING 1111 doffile not specified, using default: "<<opt.dofFile<<endl;
	}

	opt.selection = OP.getString("selection");
	if (OP.fail()){
	  opt.selection = "resn ASX+ASN+NAG+BMA+MAN";
	}

	cout << "Done setting options\n";
	return (opt);

}


bool fixAtomNames(System &_sys){
  
  bool returnValue = false;


  for (uint i = 0; i < _sys.positionSize();i++){
    if (_sys.getPosition(i).getResidueName() != "MAN" ||
	_sys.getPosition(i).getResidueName() != "BMA" ||
	_sys.getPosition(i).getResidueName() != "NAG"){
      continue;
    }

    MSLOUT.stream () << "Working on residue: "<<_sys.getPosition(i).toString()<<endl;
    // Get neighbors
    
    // Get linkages

    // Rename the oxygens
    
  }

  return (returnValue);
}
