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

#include "BackRub.h"
#include "Transforms.h"
#include "RandomNumberGenerator.h"
#include "PDBWriter.h"
#include "PyMolVisualization.h"

using namespace MSL;
using namespace std;

#include "MslOut.h"

static MslOut MSLOUT("BackRub");

BackRub::BackRub(){
	// Random number generator
	rng.setTimeBasedSeed();
	double seed = rng.getSeed();
	MSLOUT.stream() << "SEED: "<<seed<<endl;
}

BackRub::~BackRub(){
}

std::vector<AtomContainer *> BackRub::multiSample(Chain &_ch, int _startResIndex, int _endResIndex,  int _numSamplingAttempts, int _numModels, double _rmsdCutoff) {

	vector<AtomContainer *> theList;

	Chain origChain(_ch);

		
	uint m = 0;
	do {
		MSLOUT.fprintf(stdout, "NEW LOOP\n");
		Chain tmpChain(origChain);

		for (uint s = 0 ; s < _numSamplingAttempts;s++){

			// Rub the 3 residue segment modifying those atoms in tmpChain
			Chain tmpChain2(tmpChain);
			MSLOUT.fprintf(stdout, "Rubbing %d\n",s);
			bool rubOk = doRub(tmpChain2,_startResIndex,_endResIndex);
			MSLOUT.fprintf(stdout, "Rubbed\n");
			if (!rubOk){
				MSLOUT.fprintf(stdout, "Rub failed\n");
				continue;
			} 
			MSLOUT.fprintf(stdout, "Copying Chain\n");
			tmpChain = tmpChain2;
		}

		MSLOUT.fprintf(stdout, "Copying tmpChain into AC\n");

		// Create AtomContainer from tmpChain
		AtomContainer *ac = new AtomContainer(tmpChain.getAtomPointers());

		MSLOUT.fprintf(stdout, "Results push_back\n");
		// Check RMSD (not implemented yet)

		theList.push_back(ac);
		MSLOUT.fprintf(stdout, "Results pushed\n");
		m++;
	} while (m < _numModels);

	MSLOUT.fprintf(stdout, "Returning theList\n");
	return theList;
}

bool BackRub::doRub(Chain &_ch, int _startRes, int _endRes) {

		// Transforms object
		Transforms t;

		// Randomly pick a 3 residue rotation vector, within our bounds of stem1,stem2
		int randRange=(_endRes - _startRes -3 );
		int randInt  = rng.getRandomInt(randRange);
		int startRes = _startRes + randInt;
		int endRes   = startRes + 2;
		
		//CartesianPoint mainRotVector = stem2("CA").getCoor() - stem1("CA").getCoor();		
		Residue &res1 = _ch.getResidue(startRes);
		Residue &res2 = _ch.getResidue(startRes+1);
		Residue &res3 = _ch.getResidue(endRes);

		MSLOUT.fprintf(stdout, "Residues: %s %s %s\n",res1.getPositionId().c_str(),res2.getPositionId().c_str(),res3.getPositionId().c_str());

		// Weird things happen when PRO is at residue 2, have not investigated the root of the problem yet. For now just skip them.
		if (res1.getResidueName() == "PRO" || res2.getResidueName() == "PRO" || res3.getResidueName() == "PRO"){
		  return false;
		}	
	
		if (!res1.atomExists("CA")  || !res2.atomExists("CA") || !res3.atomExists("CA")){
			cerr << "ERROR BackRub::localSample() res1,res2,res3 for defining mainRotationAxis does not have CA atom.\n";
			exit(1);
		}

		if (!res2.atomExists("O")  || !res3.atomExists("O")){
			cerr << "ERROR BackRub::localSample() res2,res3 for defining mainRotationAxis does not have O atom.\n";
			exit(1);
		}


		
		/*
		  Store some atom coordinates before major rotation.  We will use carboxyl oxygens for the minor rotations.
		     That is:     O1-Calpha1-Calpha2-O2

		     This seemed to work better than: O1-Calpha1-Calpha2-N2
		 */
		CartesianPoint preO1 = res1("O").getCoor();
		CartesianPoint preO2 = res2("O").getCoor();


		// For second minor rotation we will try to use Oxygen of endRes+1 residue, if it doesn't exist use Nitrogen of endRes.
		CartesianPoint preAt2 = res2("N").getCoor();
		if (endRes+1 < _ch.positionSize() && _ch.getResidue(endRes+1).atomExists("O")){
			
			preAt2 = _ch.getResidue(endRes+1)("O").getCoor();
		}


		CartesianPoint mainRotVector = res3("CA").getCoor();
		double omega = 0.0;
		while (omega == 0.0){
			omega = rng.getRandomInt(20);
		}

		// assign  +/- randomly.
		if (rng.getRandomDouble() > 0.5){
			omega *= -1;
		}


		MSLOUT.fprintf(stdout, "\nBackRub: %s %3d %3s, %3d %3s : omega: %8.3f", res1.getChainId().c_str(), res1.getResidueNumber(), res1.getResidueName().c_str(), res3.getResidueNumber(), res3.getResidueName().c_str(),omega);
		// Rotate by omega

		// Rotate C=0, of res1 residue, by omega
		t.rotate(res1("C"),omega,mainRotVector, res1("CA").getCoor());
		t.rotate(res1("O"),omega,mainRotVector, res1("CA").getCoor());


		// Rotate all atoms of residues in between stems by omega  (for loop if larger than 3 residue fragment rotations...)
		t.rotate(res2.getAtomPointers(), omega, mainRotVector, res1("CA").getCoor());


		// Rotate N-H, of stem2 residue, by omega
		t.rotate(res3("N"),omega,mainRotVector, res1("CA").getCoor());
		if (res3.atomExists("H")){
			t.rotate(res3("H"),omega,mainRotVector, res1("CA").getCoor());
		}

		/*
		char pname[80];
		sprintf(pname,"/tmp/frag-%04d.pdb",f);
		PDBWriter pout;
		pout.open(pname);
		pout.write(_ch.getAtomPointers());
		pout.close();
		*/

		// ROTATE PEPTIDE BOND ATOMS BACK AS CLOSE TO ORIGINAL POSITION AS POSSIBLE

		// Do compensatory rotation along res1->res2 vector, try to put "O" atom back onto preO1. using O1-CA1-CA2-N1 dihedral angle
		//   this also should reduce strain placed on N-Calpha-C' bond angle due to initial rotation
		MSLOUT.fprintf(stdout, " minor1 ");
		doMinorRotation(res1,res2,preO1,preO2);

		/*
		char pname2[80];
		sprintf(pname2,"/tmp/minor1-%04d.pdb",f);
		pout.open(pname2);
		pout.write(_ch.getAtomPointers());
		pout.close();
		*/

		// Do compensatory rotation along res2->res3 vector, try to put "O" atom back onto preO2. using O2-CA2-CA3-N2 dihedral angle
		//   this also should reduce strain placed on C'-Calpha-N bond angle due to initial rotation
		MSLOUT.fprintf(stdout, " minor2 ");
		doMinorRotation(res2,res3,preO2,preAt2);
		MSLOUT.fprintf(stdout,"\n");
		/*
		char pname3[80];
		sprintf(pname3,"/tmp/minor2-%04d.pdb",f);
		pout.open(pname3);
		pout.write(_ch.getAtomPointers());
		pout.close();
		*/

		// Quick check. make sure the Ca-Ca distances are still ok.
		double dist1 = res1("CA").distance(res2("CA"));
		double dist2 = res2("CA").distance(res3("CA"));

		if (fabs(dist1 - 3.8) > 0.1 || fabs(dist2 - 3.8) > 0.1) return false;

		return true;
}

void BackRub::localSample(Chain &_ch, int _startResIndex, int _endResIndex, int _numFragments){

	Residue &stem1 = _ch.getResidue(_startResIndex);
	Residue &stem2 = _ch.getResidue(_endResIndex);

	if (!stem1.atomExists("CA") || !stem2.atomExists("CA")){
		cerr << "ERROR BackRub::localSample() one or both stem residues does not contain a 'CA' atom."<<endl;	
		exit(1);
	}

	
	RandomNumberGenerator rng;
	rng.setTimeBasedSeed();
	double seed = rng.getSeed();
	MSLOUT.stream() << "SEED: "<<seed<<endl;

	stringstream ss;
	Transforms t;
	for (uint f = 0; f < _numFragments;f++){

		// Randomly pick a 3 residue rotation vector, within our bounds of stem1,stem2
		int randRange=(_endResIndex - _startResIndex -3 );
		int randInt  = rng.getRandomInt( randRange);
		int startRes = _startResIndex + randInt;
		int endRes   = startRes + 2;
		
		//CartesianPoint mainRotVector = stem2("CA").getCoor() - stem1("CA").getCoor();		
		Residue &res1 = _ch.getResidue(startRes);
		Residue &res2 = _ch.getResidue(startRes+1);
		Residue &res3 = _ch.getResidue(endRes);

		// Weird things happen when PRO is at residue 2, have not investigated the root of the problem yet. For now just skip them.
		if (res2.getResidueName() == "PRO"){
		  continue;
		}	
	
		if (!res1.atomExists("CA")  || !res2.atomExists("CA") || !res3.atomExists("CA")){
			cerr << "ERROR BackRub::localSample() res1,res2,res3 for defining mainRotationAxis does not have CA atom.\n";
			exit(1);
		}

		if (!res2.atomExists("O")  || !res3.atomExists("O")){
			cerr << "ERROR BackRub::localSample() res2,res3 for defining mainRotationAxis does not have O atom.\n";
			exit(1);
		}


		
		/*
		  Store some atom coordinates before major rotation.  We will use carboxyl oxygens for the minor rotations.
		     That is:     O1-Calpha1-Calpha2-O2

		     This seemed to work better than: O1-Calpha1-Calpha2-N2
		 */
		CartesianPoint preO1 = res1("O").getCoor();
		CartesianPoint preO2 = res2("O").getCoor();


		// For second minor rotation we will try to use Oxygen of endRes+1 residue, if it doesn't exist use Nitrogen of endRes.
		CartesianPoint preAt2 = res2("N").getCoor();
		if (endRes+1 < _ch.positionSize() && _ch.getResidue(endRes+1).atomExists("O")){
			
			preAt2 = _ch.getResidue(endRes+1)("O").getCoor();
		}


		CartesianPoint mainRotVector = res3("CA").getCoor();
		double omega = 0.0;
		while (omega == 0.0){
			omega = rng.getRandomInt(20);
		}

		// assign  +/- randomly.
		if (rng.getRandomDouble() > 0.5){
			omega *= -1;
		}


		MSLOUT.fprintf(stdout, "\nBackRub%04d: %s %3d %3s, %3d %3s : omega: %8.3f",f, res1.getChainId().c_str(), res1.getResidueNumber(), res1.getResidueName().c_str(), res3.getResidueNumber(), res3.getResidueName().c_str(),omega);
		// Rotate by omega

		// Rotate C=0, of res1 residue, by omega
		t.rotate(res1("C"),omega,mainRotVector, res1("CA").getCoor());
		t.rotate(res1("O"),omega,mainRotVector, res1("CA").getCoor());


		// Rotate all atoms of residues in between stems by omega  (for loop if larger than 3 residue fragment rotations...)
		t.rotate(res2.getAtomPointers(), omega, mainRotVector, res1("CA").getCoor());


		// Rotate N-H, of stem2 residue, by omega
		t.rotate(res3("N"),omega,mainRotVector, res1("CA").getCoor());
		if (res3.atomExists("H")){
			t.rotate(res3("H"),omega,mainRotVector, res1("CA").getCoor());
		}

		/*
		char pname[80];
		sprintf(pname,"/tmp/frag-%04d.pdb",f);
		PDBWriter pout;
		pout.open(pname);
		pout.write(_ch.getAtomPointers());
		pout.close();
		*/

		// ROTATE PEPTIDE BOND ATOMS BACK AS CLOSE TO ORIGINAL POSITION AS POSSIBLE

		// Do compensatory rotation along res1->res2 vector, try to put "O" atom back onto preO1. using O1-CA1-CA2-N1 dihedral angle
		//   this also should reduce strain placed on N-Calpha-C' bond angle due to initial rotation
		MSLOUT.fprintf(stdout, " minor1 ");
		doMinorRotation(res1,res2,preO1,preO2);

		/*
		char pname2[80];
		sprintf(pname2,"/tmp/minor1-%04d.pdb",f);
		pout.open(pname2);
		pout.write(_ch.getAtomPointers());
		pout.close();
		*/

		// Do compensatory rotation along res2->res3 vector, try to put "O" atom back onto preO2. using O2-CA2-CA3-N2 dihedral angle
		//   this also should reduce strain placed on C'-Calpha-N bond angle due to initial rotation
		MSLOUT.fprintf(stdout, " minor2 ");
		doMinorRotation(res2,res3,preO2,preAt2);
		MSLOUT.fprintf(stdout,"\n");
		/*
		char pname3[80];
		sprintf(pname3,"/tmp/minor2-%04d.pdb",f);
		pout.open(pname3);
		pout.write(_ch.getAtomPointers());
		pout.close();
		*/

		// Quick check. make sure the Ca-Ca distances are still ok.
		double dist1 = res1("CA").distance(res2("CA"));
		double dist2 = res2("CA").distance(res3("CA"));

		if (fabs(dist1 - 3.8) > 0.1 || fabs(dist2 - 3.8) > 0.1) continue;


		sys.addAtoms(_ch.getAtomPointers());

		// Write out to stringstream at this point...
		ss << "MODEL"<<endl;
		PDBWriter pss;
		pss.open(ss);
		pss.write(_ch.getAtomPointers());
		pss.close();		
		ss << "ENDMDL"<<endl;

	}

	sysNMRFormat = ss.str();

}

void BackRub::doMinorRotation(Residue &_r1, Residue &_r2, CartesianPoint &_targetAtom1, CartesianPoint &_targetAtom2){



	// Define minor rotation vector	
	CartesianPoint minorRotVector = _r2("CA").getCoor();

	// Get pre-dihedral
	double preDihedral = _targetAtom1.dihedral(_r1("CA").getCoor(),_r2("CA").getCoor(),_targetAtom2);

	// Get post-dihedral
	double postDihedral = _r1("O").getCoor().dihedral(_r1("CA").getCoor(),_r2("CA").getCoor(),_targetAtom2);

	// Rotate back by post-pre
	double alpha = postDihedral - preDihedral;
	MSLOUT.fprintf(stdout, "preDihedral = %8.3f postDihedral = %8.3f alpha = %8.3f", preDihedral, postDihedral, alpha);




	// Create a local transform object to rotate atoms...
	Transforms t;
	t.rotate(_r1("O"),alpha, minorRotVector, _r1("CA").getCoor());
	t.rotate(_r1("C"),alpha, minorRotVector, _r1("CA").getCoor());
	t.rotate(_r2("N"),alpha, minorRotVector, _r1("CA").getCoor());
	if (_r2.atomExists("H")){
		t.rotate(_r2("H"),alpha, minorRotVector, _r1("CA").getCoor());
	}


}


