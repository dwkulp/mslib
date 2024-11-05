#include <iostream>
#include <cstdlib>
#include <queue>
#include "signal.h"
#include "System.h"
#include "Frame.h"
#include "Timer.h"
#include "RegEx.h"
#include "Transforms.h"
#include "PDBWriter.h"
#include "SasaCalculator.h"
#include "AtomContainer.h"
#include "AtomPointerVector.h"
#include "OptionParser.h"
#include "MslOut.h"
#include "MslExceptions.h"


#include "searchForNanos.h"



using namespace std;
using namespace MSL;


// MslOut 
static MslOut MSLOUT("searchForNanos");


int main(int argc, char *argv[]) {
  
  // Parse commandline options
  Options opt = setupOptions(argc,argv);

  // MslOut can suppress output, to make example output clean
  //MSLOUT.turnOn("searchForNanos");

  vector<string> lines;
  MslTools::readTextFile(lines,opt.pdbs);
  

  for (uint i = 0; i < lines.size();i++){

    System sys;
    sys.readStructureFile(lines[i]);

    if (sys.chainSize() < opt.numChains) continue;

    AtomContainer ac;
    for (uint c = 0; c < sys.chainSize(); c++){
      if (sys.getChain(c).getPosition(0).atomExists("CA")){
	//fprintf(stdout,"Adding: %s\n",sys.getChain(c).getPosition(0).toString().c_str());
	ac.addAtom(sys.getChain(c).getPosition(0).getAtom("CA"));
	//fprintf(stdout,"Done\n");
      }
    }

    vector<std::pair<int,int> > closeAtoms;
    map<int,int> closeFlags;
    for (uint ii = 0; ii < ac.size();ii++){
      for (uint jj = ii+1;jj < ac.size();jj++){
	//fprintf(stdout,"ii jj dist: %s %s %6.2f\n",ac(ii).getAtomId().c_str(),ac(jj).getAtomId().c_str(),ac(ii).distance(ac(jj)));
	if (ac(ii).distance(ac(jj)) < opt.distCutoff){
	  closeAtoms.push_back(std::pair<int,int>(ii,jj));
	  if (closeFlags.find(ii) == closeFlags.end()){
	    closeFlags[ii] = 1;
	  } else {
	    closeFlags[ii] = closeFlags[ii]+1;
	  }
	  if (closeFlags.find(jj) == closeFlags.end()){
	    closeFlags[jj] = 1;
	  } else {
	    closeFlags[jj] = closeFlags[jj]+1;
	  }
	  //fprintf(stdout, "close: %d %d %d %d\n", ii,jj,closeFlags[ii],closeFlags[jj]);
	}
      }
    }

    

    // Output any that meet requirment
    for (auto it = closeFlags.cbegin(); it != closeFlags.cend();++it){
      
      if (it->second == opt.numCloseAts-1){
	fprintf(stdout, "DATA NTERM %s %d %s\n", lines[i].c_str(),(uint)closeAtoms.size(),ac(it->first).getAtomId().c_str());
      }
    }


    ac.removeAllAtoms();

    for (uint c = 0; c < sys.chainSize(); c++){
      if (sys.getChain(c).getPosition(sys.getChain(c).positionSize()-1).atomExists("CA")){
	//fprintf(stdout,"Adding: %s\n",sys.getChain(c).getPosition(sys.getChain(c).positionSize()-1).toString().c_str());
	ac.addAtom(sys.getChain(c).getPosition(sys.getChain(c).positionSize()-1).getAtom("CA"));
	//fprintf(stdout,"Done\n");
      }
    }
    
    closeAtoms.clear();
    closeFlags.clear();

    for (uint ii = 0; ii < ac.size();ii++){
      for (uint jj = ii+1;jj < ac.size();jj++){
	//fprintf(stdout,"ii jj dist: %s %s %6.2f\n",ac(ii).getAtomId().c_str(),ac(jj).getAtomId().c_str(),ac(ii).distance(ac(jj)));
	if (ac(ii).distance(ac(jj)) < opt.distCutoff){
	  closeAtoms.push_back(std::pair<int,int>(ii,jj));
	  if (closeFlags.find(ii) == closeFlags.end()){
	    closeFlags[ii] = 1;
	  } else {
	    closeFlags[ii] = closeFlags[ii]+1;
	  }
	  if (closeFlags.find(jj) == closeFlags.end()){
	    closeFlags[jj] = 1;
	  } else {
	    closeFlags[jj] = closeFlags[jj]+1;
	  }
	  //fprintf(stdout, "close: %d %d %d %d\n", ii,jj,closeFlags[ii],closeFlags[jj]);
	}
      }
    }

    for (auto it = closeFlags.cbegin(); it != closeFlags.cend();++it){
      
      if (it->second == opt.numCloseAts-1){
	fprintf(stdout, "DATA CTERM %s %d %s\n", lines[i].c_str(),(uint)closeAtoms.size(),ac(it->first).getAtomId().c_str());
      }
    }

  }

	

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
    cout << "searchForNanos --pdbs PDB_LIST\n";
    exit(0);
  }

  opt.pdbs = OP.getString("pdbs");
  if (OP.fail()){
    cerr << "ERROR 1111 pdbs not specified.\n";
    exit(1111);
  }

  opt.numChains = OP.getInt("numChains");
  if (OP.fail()){
    opt.numChains = 12;
    cerr << "WARNING 1111 numChains defaulted to "<<opt.numChains<<endl;
  }

  opt.numCloseAts = OP.getInt("numCloseAts");
  if (OP.fail()){
    opt.numCloseAts = 4;
    cerr << "WARNING 1111 numCloseAts defaulted to "<<opt.numCloseAts<<endl;
  }

  opt.distCutoff = OP.getDouble("distCutoff");
  if (OP.fail()){
    opt.distCutoff = 10;
    cerr << "WARNING 1111 distCutoff defaulted to "<<opt.distCutoff<<endl;
  }

  return opt;
}

