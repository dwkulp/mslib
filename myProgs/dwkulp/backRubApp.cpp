#include <iostream>
#include <cstdlib>
#include <queue>
#include <algorithm>
#include <map>

#include "System.h"
#include "OptionParser.h"
#include "MslOut.h"
#include "SysEnv.h"
#include "SasaCalculator.h"
#include "PhiPsiStatistics.h"
#include "MonteCarloManager.h"
#include "Quench.h"
#include "PDBTopology.h"
#include "AtomSelection.h"
#include "PyMolVisualization.h"
#include "BackRub.h"

#include "backRubApp.h"

using namespace std;
using namespace MSL;

static MslOut MSLOUT("backrubApp");
static SysEnv SYSENV;

int main(int argc, char *argv[]) {

  // Parse commandline options
  Options opt = setupOptions(argc,argv);

  // Read PDB structure
  System sys;
  sys.readPdb(opt.pdb);

  BackRub br;
  vector<AtomContainer *> results = br.multiSample(sys.getChain(0), 1, sys.getChain(0).positionSize()-1, opt.numSampling, opt.numModels);

  for (uint i = 0; i < results.size(); i++){
    char name[80];
    sprintf(name,"%s_BR%06d.pdb",MslTools::getFileName(opt.pdb).c_str(),i+1);
    MSLOUT.stream() << "WRITING: "<<name<<endl;
    results[i]->writePdb((string)name);
  }

  /*
  br.localSample(sys.getChain(0),1,sys.getChain(0).positionSize()-1,opt.numModels);




  int numConfs = br.getAtomPointers().getMaxAltConf();
  AtomContainer ats(br.getAtomPointers());
  string fname = MslTools::getFileName(opt.pdb);
  for (uint i = 0; i < numConfs;i++){
    ats.setActiveConformation(i);
    char name[80];
    sprintf(name,"%s_BR%06d.pdb",fname.c_str(),i);
    MSLOUT.stream() << "WRITING: "<<name<<endl;
    ats.writePdb((string)name);
  }
  */

}




Options setupOptions(int theArgc, char * theArgv[]){
  Options opt;

  OptionParser OP;


  OP.setRequired(opt.required);
  OP.setAllowed(opt.optional);
  OP.setDefaultArguments(opt.defaultArgs); // a pdb file value can be given as a default argument without the --pdbfile option
  OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"
  OP.readArgv(theArgc, theArgv);

  if (OP.countOptions() == 0){
    cout << "Usage:" << endl;
    cout << endl;
    cout << "resurfaceSaltBridges --pdb PDB\n";

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
    cerr << "ERROR 1111 pdb not specified.\n";
    exit(1111);
  }

  opt.numModels = OP.getInt("numModels");
  if (OP.fail()){
    opt.numModels = 100;
  }

  opt.numSampling = OP.getInt("numSampling");
  if (OP.fail()){
    opt.numSampling = 10;
  }

  MSLOUT.stream() << "Options:\n"<<OP<<endl;
  return opt;
}
