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


#include "transformSelections.h"

using namespace std;
using namespace MSL;


// MslOut 
static MslOut MSLOUT("transformSelections");
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
  AtomPointerVector ats3    = sel.select(opt.sele3);	
  AtomPointerVector ats_all = ats1 + ats2 + ats3;
  
  AtomSelection sel1(ats1);
  AtomPointerVector ats1_ca = sel1.select("name CA");
  
  AtomSelection sel2(ats2);
  AtomPointerVector ats2_ca = sel2.select("name CA");

  AtomSelection sel3(ats3);
  AtomPointerVector ats3_ca = sel3.select("name CA");

  AtomContainer acopy(sys.getAtomPointers());
  
  Transforms tm;
  tm.rmsdAlignment(ats1_ca,ats3_ca, acopy.getAtomPointers());
  acopy.writePdb("out.pdb");

  /*  
  // Align 1 to 3
  // Apply align orig 1 to 2
  // Align 2 to 3
  // Apply align orig 2 to 3
  tm.rmsdAlign(ats1,ats3);
  tm.rmsdAlign(ats1_tm,ats2_tm,ats1);
  tm.rmsdAlign(ats2_tm,ats3_tm,ats1);
  */
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
    cout << "transformSelections --pdb pdb --sele1 sele1 --sele2 sele2 --sele3 sele3 --numRepeats 1\n";

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

  opt.sele3 = OP.getString("sele3");
  if (OP.fail()){
    cerr << "ERROR 1111 sele3 not specified.\n";
    exit(1111);
  }
  
  opt.numRepeats = OP.getInt("numRepeats");
  if (OP.fail()){
    cerr << "ERROR 1111 numRepeats not specified.\n";
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
