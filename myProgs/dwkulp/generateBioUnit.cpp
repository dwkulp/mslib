#include <iostream>
#include <cstdlib>
#include <queue>

#include "System.h"
#include "Transforms.h"
#include "PDBWriter.h"
#include "AtomPointerVector.h"
#include "OptionParser.h"
#include "MslOut.h"

#include "generateBioUnit.h"


using namespace std;
using namespace MSL;


// MslOut 
static MslOut MSLOUT("generateBioUnit");



int main(int argc, char *argv[]) {

  // Parse commandline options
  Options opt = setupOptions(argc,argv);

  // Read in the reference pdb file
  System sys;
  sys.readPdb(opt.pdb);

  // Get Biological Unit Operations
  vector<Matrix *> bioUnitRotations            = sys.getPDBReader()->getBioUnitRotations();
  vector<CartesianPoint *> bioUnitTranslations = sys.getPDBReader()->getBioUnitTranslations();

  // Transforms object 
  Transforms tr;

  // Setup a PDBWriter
  PDBWriter pout;
  pout.open(opt.outPdb);
  pout.write(sys.getAtomPointers(),false,false,true); // write input coordinates as first model, do I need to do this or is this unit encoded in the BIOMT?

  // For each matrix (or model/unit)
  for (uint m = 0; m < bioUnitRotations.size();m++){
    MSLOUT.stream()<< "Model "<<m<<endl;

    // Rotate
    tr.rotate(sys.getAtomPointers(),*bioUnitRotations[m]);

    // Translate
    tr.translate(sys.getAtomPointers(),*bioUnitTranslations[m]);

    // Write model
    pout.write(sys.getAtomPointers(),false,false,true);
    
  }
  pout.close();
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
    cout << "generateBioUnit --pdb PDB\n";
    exit(0);
  }


  opt.configfile = OP.getString("configfile");
  if (opt.configfile != "") {
    OP.readFile(opt.configfile);
    if (OP.fail()) {
      cerr << "ERROR couldn't read : "<<opt.configfile<<endl;
      exit(1);
    }
  }


  opt.pdb = OP.getString("pdb");
  if (OP.fail()){
    cerr << "ERROR 1111 pdb not specified.\n";
    exit(1111);
  }


  opt.outPdb = OP.getString("outPdb");
  if (OP.fail()){
    opt.outPdb = "foo.pdb";
    cerr << "WARNING --outPdb is set to : "<<opt.outPdb<<endl;
  }

  return opt;
}
