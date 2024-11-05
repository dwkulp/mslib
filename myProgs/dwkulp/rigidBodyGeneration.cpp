#include <iostream>
#include <cstdlib>
#include <queue>
#include <algorithm>
#include <map>

#include "System.h"
#include "OptionParser.h"
#include "Frame.h"
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


#include "rigidBodyGeneration.h"

using namespace std;
using namespace MSL;


// MslOut 
static MslOut MSLOUT("rigidBodyGeneration");
static SysEnv SYSENV;


int main(int argc, char *argv[]) {

  // Parse commandline options
  Options opt = setupOptions(argc,argv);

  // Read in original pdb
  System sys;
  sys.readStructureFile(opt.pdb);
  sys.getAtomPointers().saveCoor("orig");

  AtomSelection sel(sys.getAtomPointers());
  AtomPointerVector ats1    = sel.select(opt.sele1);
  AtomPointerVector ats2    = sel.select(opt.sele2);	
  AtomPointerVector ats_all = ats1 + ats2;

  AtomSelection sel1(ats1);
  AtomPointerVector ats1_ca = sel1.select("name CA");
  
  AtomSelection sel2(ats2);
  AtomPointerVector ats2_ca = sel2.select("name CA");
  
  CartesianPoint center1    = ats1.getGeometricCenter();
  CartesianPoint center2    = ats2.getGeometricCenter();
  CartesianPoint center_all = ats_all.getGeometricCenter();

  PyMolVisualization pymol;
  pymol.createAtom(center1, "Center1");
  pymol.createAtom(center2, "Center2");
  pymol.createAtom(center_all, "CenterAll");

  Line line1(center_all,center1-center_all);
  line1.setName("line1");
  line1.setOutputFormat("pymol");
  
  Line line2(center_all,center2-center_all);
  line2.setName("line2");
  line2.setOutputFormat("pymol");

  // Now do some rigid body movements!
  
  /*
  Transforms tm;
  translate(ats1, 20, center_all, center1, line1);

  PDBWriter pout;
  pout.open(MslTools::stringf("test20-%02d.pdb",1));
  pout.write(ats1);
  pout.close();
  */
  
  Transforms tm;
  // Translate sele1 and sele2 first
  // Rotate around line1 and line2 second
  // Rotate around norm last
  int solution = 1;
  for (uint t1 = 10; t1 <= 25; t1+=5){
    
    double factor1 = (t1-(center1.distance(center_all)))/(center1.distance(center_all));
    CartesianPoint dir1 = line1.getCenter() - line1.getDirection()*factor1;
    tm.translate(ats1, CartesianPoint(line1.getCenter()-dir1));

    sys.getAtomPointers().saveCoor("preR1");
      for (uint r1 = 0; r1< 180;r1+=36){
	CartesianPoint axis1(line1.getDirection());
	axis1.getUnit();
	tm.rotate(ats1, r1, center_all-axis1,center_all);

	sys.getAtomPointers().saveCoor("preR2");
	for (uint r2 = 0; r2< 180;r2+=36){
	  CartesianPoint axis2(line2.getDirection());
	  axis2.getUnit();
	  tm.rotate(ats2, r2, center_all-axis2,center_all);

	  sys.getAtomPointers().saveCoor("preN");
	  for (int n = 0; n <= 120;n+=30){
	    CartesianPoint norm(line1.getDirection().cross(line2.getDirection()));
	    norm = norm.getUnit();
	    tm.rotate(ats1, n, center_all-norm,center_all);


	    // Check for clashes (more than 3 at 6 Angstroms), skip
	    int clashes = ats1_ca.clashCheck(ats2_ca,36) ;
	    if (clashes > 3){
	      continue;
	    }
	    cout << "Solution : "<<solution<< " [ "<<t1<<","<<r1<<","<<r2<<","<<n<<" ] \n";
	    // Write PDB
	    PDBWriter pout;
	    pout.open(MslTools::stringf("rigid_body_%05d.pdb",solution));
	    pout.write(ats1);
	    pout.write(ats2);
	    pout.close();

	    // apply save coor orig
	    sys.getAtomPointers().applySavedCoor("preN");

	    //pout.open(MslTools::stringf("preN_%05d.pdb",solution));
	    //pout.write(ats1);
	    //pout.write(ats2);
	    //pout.close();
	    
	    solution += 1; 
	    //if (solution == 10){
	    //  exit(1111);
	    //}
	    
	  } // for n
	  
	  sys.getAtomPointers().applySavedCoor("preR2");
	} // for r2
	sys.getAtomPointers().applySavedCoor("preR1");
      } // for r1
      sys.getAtomPointers().applySavedCoor("orig");
  } // for t1
  
  /*
  // Rotate around normal axis
  CartesianPoint norm(line1.getDirection().cross(line2.getDirection()));
  norm = norm.getUnit();
  tm.rotate(ats1, 45, center_all-norm,center_all);

  PDBWriter pout;
  pout.open(MslTools::stringf("testNorm-%02d.pdb",1));
  pout.write(ats1);
  pout.close();

  pymol.createAtom(ats1.getGeometricCenter(),"GC1_norm"); 
  */  
  /*
  // Rotate along a line
  CartesianPoint axis1(line1.getDirection());
  axis1.getUnit();
  tm.rotate(ats1, 90, center_all-axis1,center_all);

  PDBWriter pout;
  pout.open(MslTools::stringf("testRotate-%02d.pdb",1));
  pout.write(ats1);
  pout.close();
  */
  
  /*
  // Change distances by translating along a line  (how do I covert this factor into a distance)
  double dist_sele1 = 10;
  double factor = (dist_sele1+(center1.distance(center_all)))/(center1.distance(center_all));
  double factor2 = (dist_sele1-(center1.distance(center_all)))/(center1.distance(center_all));
  
  CartesianPoint dir1a = line1.getCenter() + line1.getDirection()*1.75;
  CartesianPoint dir1b = line1.getCenter() - line1.getDirection()*1.75;
  CartesianPoint dir1c = line1.getCenter() + line1.getDirection()*factor;
  CartesianPoint dir1d = line1.getCenter() - line1.getDirection()*factor2;
  
  pymol.createAtom(dir1a,"dir1a"); // is the opposite direction along line than I want to translate.
  pymol.createAtom(dir1c,"dir1c"); // is the opposite direction along line than I want to translate.
  
  ats1.saveCoor("orig");
  tm.translate(ats1, CartesianPoint(line1.getCenter()-dir1b));
  PDBWriter pout;
  pout.open(MslTools::stringf("testB-%02d.pdb",1));
  pout.write(ats1);
  pout.close();
  
  ats1.applySavedCoor("orig");
  tm.translate(ats1, CartesianPoint(line1.getCenter()-dir1d));
  pout.open(MslTools::stringf("testD-%02d.pdb",1));
  pout.write(ats1);
  pout.close();

  pymol.createAtom(ats1.getGeometricCenter(),"GC1d"); 
  */
  
  ofstream fout;
  fout.open(MslTools::stringf("stuff.py"));
  fout << pymol;
  fout << line1.toString();
  fout << line2.toString();
  fout.close();
    
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
    cout << "rigidBodyGeneration --pdb pdb --sele1 sele1 --sele2 sele2\n";

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

  opt.sele1 = OP.getString("sele1");
  if (OP.fail()){
    cerr << "ERROR 1111 sele1 not specified.\n";
    exit(1111);
  }
  
  opt.sele2 = OP.getString("sele2");
  if (OP.fail()){
    cerr << "ERROR 1111 sele2 not specified.\n";
    exit(1111);
  }
  
  MSLOUT.stream() << "Options:\n"<<OP<<endl;
  return opt;
}


void translate(AtomPointerVector &_ats, double _dist, CartesianPoint _center, CartesianPoint _GC, Line _line){
  double factor = (_dist-(_GC.distance(_center)))/(_GC.distance(_center));
  CartesianPoint dir = _line.getCenter() - _line.getDirection()*factor;

  Transforms tm;
  tm.translate(_ats, CartesianPoint(_line.getCenter()-dir));
  
}
