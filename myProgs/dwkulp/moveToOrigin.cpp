#include <iostream>
#include <cstdlib>
#include <queue>
#include <algorithm>
#include <map>

#include "System.h"
#include "OptionParser.h"
#include "Frame.h"
#include "Line.h"
#include "MslOut.h"
#include "SysEnv.h"
#include "SasaCalculator.h"
#include "PhiPsiStatistics.h"
#include "MonteCarloManager.h"
#include "Quench.h"
#include "PDBTopology.h"
#include "AtomSelection.h"
#include "PyMolVisualization.h"
#include "AtomSelection.h"
#include "AtomContainer.h"
#include "RandomNumberGenerator.h"
#include "PDBWriter.h"
#include "Helanal.h"
#include "Transforms.h"

#include "moveToOrigin.h"

using namespace std;
using namespace MSL;


// MslOut 
static MslOut MSLOUT("moveToOrigin");
static SysEnv SYSENV;

  
int main(int argc, char *argv[]) {

  // Parse commandline options
  Options opt = setupOptions(argc,argv);

  // Read in original pdb
  System sys;
  sys.readStructureFile(opt.pdb);
  sys.getAtomPointers().saveCoor("orig");

  // Create a pymol object and global frame (should be in some function!)
  PyMolVisualization pymol;

  // Create a global frame and add to PyMol script
  Frame globalFrame;
  CartesianPoint origin((Real)0.0f, (Real)0.0f, (Real)0.0f);
  CartesianPoint xAxis((Real)1.0f, (Real)0.0f, (Real)0.0f);
  CartesianPoint yAxis((Real)0.0f, (Real)1.0f, (Real)0.0f);
  CartesianPoint zAxis((Real)0.0f, (Real)0.0f, (Real)1.0f);
  Line zLine(origin,zAxis);
  Line xLine(origin,xAxis);
  globalFrame.computeFrameFrom2Lines(zLine,xLine);
  CartesianPoint xcen = globalFrame["X"].getCenter();
  CartesianPoint xdir = globalFrame["X"].getDirection();
  CartesianPoint ycen = globalFrame["Y"].getCenter();
  CartesianPoint ydir = globalFrame["Y"].getDirection();
  CartesianPoint zcen = globalFrame["Z"].getCenter();
  CartesianPoint zdir = globalFrame["Z"].getDirection();
  pymol.createArrow(xcen,xdir,"XaxisGlobal");
  pymol.createArrow(ycen,ydir,"YaxisGlobal");
  pymol.createArrow(zcen,zdir,"ZaxisGlobal");
  CartesianPoint center(0,0,0);
  pymol.createAtom(center, "Center");

  // Use Geometric Centers
  AtomPointerVector chA = sys.getChain("A").getAtomPointers();
  AtomPointerVector chC = sys.getChain("C").getAtomPointers();
  AtomPointerVector chE = sys.getChain("E").getAtomPointers();
  //Frame gcFrame;
  //gcFrame.computeFrameFrom3Atoms(chA.getGeometricCenter(),chC.getGeometricCenter(),chE.getGeometricCenter());
  //gcFrame.transformToGlobalBasis(sys.getAtomPointers());
  
  
    
  AtomSelection sel(sys.getAtomPointers());
  int bestPos = -1;
  double minRMSD = 1000;
  Frame bestFrame;
  for (uint r = 0; r<sys.getChain("A").positionSize();r++){
    sys.applySavedCoor("orig");
    cout << "Working on "<<r<<" "<<sys.getChain("A").getPosition(r)<<endl;
    char c[100];
    sprintf(c, "resi %d and name CA", sys.getChain("A").getPosition(r).getResidueNumber());
    AtomPointerVector ats = sel.select((std::string)c);
    
    Frame aFrame;
    aFrame.computeFrameFrom3Atoms(ats(0),ats(1),ats(2),true);
    aFrame.transformToGlobalBasis(sys.getAtomPointers());


    AtomContainer chA_2(chA);
    AtomContainer chA_3(chA);
  
    Transforms tm;
    tm.rotate(chA_2.getAtomPointers(),120,globalFrame["Z"].getCenter()-globalFrame["Z"].getDirection(),globalFrame["Z"].getCenter());
    tm.rotate(chA_3.getAtomPointers(),240,globalFrame["Z"].getCenter()-globalFrame["Z"].getDirection(),globalFrame["Z"].getCenter());

    // Compute RMSD
    double rmsd1 = chA_2.getAtomPointers().rmsd(chE);
    double rmsd2 = chA_3.getAtomPointers().rmsd(chC);
    double rmsd = (rmsd1+rmsd2)/2;

    if (rmsd <= minRMSD){
      cout << "\t best yet: "<<r<<" rmsd: "<<rmsd<<endl;
      bestPos = r;
      minRMSD = rmsd;
      bestFrame = aFrame;
      sys.getAtomPointers().saveCoor("best");
    } else {
      cout << "\t not best: "<<r<<" rmsd: "<<rmsd<<endl;
    }
    
  }

  // revert to best
  sys.applySavedCoor("best");

  // write pdb
  sys.writePdb("best_orig.pdb");
  
  // create rotated chC, chE
  AtomContainer chA_1(chA);
  AtomContainer chA_2(chA);
  AtomContainer chA_3(chA);
  
  Transforms tm;
  tm.rotate(chA_2.getAtomPointers(),120,globalFrame["Z"].getCenter()-globalFrame["Z"].getDirection(),globalFrame["Z"].getCenter());
  tm.rotate(chA_3.getAtomPointers(),240,globalFrame["Z"].getCenter()-globalFrame["Z"].getDirection(),globalFrame["Z"].getCenter());

  chA_1.writePdb("chA.pdb");
  chA_2.writePdb("chE.pdb");
  chA_3.writePdb("chC.pdb");
  
  // write aFrame to pymol
  xcen = bestFrame["X"].getCenter();
  xdir = bestFrame["X"].getDirection();
  ycen = bestFrame["Y"].getCenter();
  ydir = bestFrame["Y"].getDirection();
  zcen = bestFrame["Z"].getCenter();
  zdir = bestFrame["Z"].getDirection();
  pymol.createArrow(xcen,xdir,"XaxisFrame");
  pymol.createArrow(ycen,ydir,"YaxisFrame");
  pymol.createArrow(zcen,zdir,"ZaxisFrame");

  ofstream fout;
  fout.open(MslTools::stringf("axisGlobalOriginBest.py"));
  fout << pymol;
  fout.close();
  
  //double distA = chA.getGeometricCenter().distance(center);
  //double distC = chC.getGeometricCenter().distance(center);
  //double distE = chE.getGeometricCenter().distance(center);
  //
  //MSLOUT.stream() << "A: "<<distA<<" C: "<<distC<<" E: "<<distE<<endl;


  MSLOUT.stream() << "Done."<<endl;
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
    cout << "simpleNanoparticleModeling --pdb pdb\n";

    cout << "\nprogram options: "<<endl;
    for (uint i = 0; i < opt.required.size();i++){
      cout <<"R  --"<<opt.required[i]<<"  "<<endl;
    }
    cout <<endl;
    for (uint i = 0; i < opt.optional.size();i++){
      cout <<"O  --"<<opt.optional[i]<<"  "<<endl;
    }
    cout << endl;
    exit(0);
  }

  opt.pdb = OP.getString("pdb");
  if (OP.fail()){
    cerr << "ERROR 1111 pdb/cif not specified.\n";
    exit(1111);
  }

  return opt;
}



