#include <iostream>
#include <cstdlib>

#include "AtomContainer.h"
#include "AtomPointerVector.h"
#include "MslTools.h"
#include "PDBTopology.h"
#include "PDBWriter.h"
#include "System.h"
#include "MslOut.h"
#include "PyMolVisualization.h"
#include "AtomSelection.h"
//#include "VectorHashing.h"
#include "OptionParser.h"
#include "Transforms.h"
#include "Timer.h"
#include "VectorHashing.h"
#include "Frame.h"

#include "alignByFrame.h"

using namespace std;
using namespace MSL;

// MslOut 
static MslOut MSLOUT("alignByFrame");

int main(int argc, char *argv[]) {
        Options opt = setupOptions(argc,argv);

	AtomContainer ref;
	ref.readPdb(opt.pdb1);

	AtomSelection refSel(ref.getAtomPointers());
        AtomPointerVector refSelAts = refSel.select(opt.sel1);
        Frame refFrame;
        if (refSelAts.size() == 3) {
          MSLOUT.stream()<<" Ref Frame from 3 Atoms\n";  
          CartesianPoint p1 = refSelAts(0).getCoor();
          CartesianPoint p2 = refSelAts(1).getCoor();
          CartesianPoint p3 = refSelAts(2).getCoor();
          refFrame.computeFrameFrom3Points(p1,p2,p3,true);
        } else {
           refFrame.computeFrameFromPCA(refSelAts);
       }

	AtomContainer pdb;
	pdb.readPdb(opt.pdb2);
		
	AtomSelection pdbSel(pdb.getAtomPointers());
        AtomPointerVector pdbSelAts = pdbSel.select(opt.sel2);

	Frame pdbFrame;
        if (pdbSelAts.size() == 3) {
          MSLOUT.stream()<<" Ref Frame from 3 Atoms\n";  
          CartesianPoint p1 = refSelAts(0).getCoor();
          CartesianPoint p2 = refSelAts(1).getCoor();
          CartesianPoint p3 = refSelAts(2).getCoor();
          pdbFrame.computeFrameFrom3Points(p1,p2,p3,true);
        } else {
          pdbFrame.computeFrameFromPCA(refSelAts);
       }

        pdbFrame.transformToGlobalBasis(pdb.getAtomPointers());
        //pdbFrame.transformAtoms(pdb.getAtomPointers(), refFrame, pdbFrame);
        refFrame.transformFromGlobalBasis(pdb.getAtomPointers());

	pdb.writePdb("out.pdb");

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
		cout << "alignByFrame --pdb1 PDB --sel1 resi 1-3 --pdb2 --sel2 resi 1-3 \n";
		exit(0);
	}
	opt.pdb1 = OP.getString("pdb1");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}

	opt.sel1 = OP.getString("sel1");
	if (OP.fail()){
		cerr << "ERROR 1111 sel not specified.\n";
		exit(1111);
	}

	opt.pdb2 = OP.getString("pdb2");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}

	opt.sel2 = OP.getString("sel2");
	if (OP.fail()){
		cerr << "ERROR 1111 sel not specified.\n";
		exit(1111);
	}

	return opt;
}
