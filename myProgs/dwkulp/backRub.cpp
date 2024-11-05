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

#include "getSSE.h"

using namespace std;
using namespace MSL;

static MslOut MSLOUT("getSSE");
static SysEnv SYSENV;

int main(int argc, char *argv[]) {

  // Parse commandline options
  Options opt = setupOptions(argc,argv);

  // Read PDB structure
  System sys;
  sys.readPdb(opt.pdb);

  map<string,int> sse_map;
  stringstream allSSE;
  stringstream previousSSE;
  string lastSSE = "";
  for (uint p = 1; p < sys.positionSize()-1;p++){
    if (!sys.getPosition(p).atomExists("CA")) continue;
    //if (!sys.getPosition(p).getAtom("CA").getSelectionFlag("resurface")) continue;

    double phi = PhiPsiStatistics::getPhi(sys.getPosition(p-1).getCurrentIdentity(),sys.getPosition(p).getCurrentIdentity());
    double psi = PhiPsiStatistics::getPsi(sys.getPosition(p).getCurrentIdentity(),sys.getPosition(p+1).getCurrentIdentity());


    string sse = "";
      if (phi < -35 &&
	  phi > -90 &&
	  psi > -70 &&
	  psi < 0){
	sse = "H";
      } else  {
	/*
	   PHI/PSI Defintion of Beta.....
	   if (phi < -95 &&
	   phi > -155 &&
	   psi > 95 &&
	   psi < 155){
	*/

	// Hbond pattern defintion
	bool isSheet = checkBetaSheet(sys, p);
	if (isSheet){
	  sse = "E";
	} else {
	  sse = "O";
	}
      }
      // SSE Boundary
      if (sse != lastSSE){
	if (previousSSE.str().length() > 4){
	  sse_map[lastSSE]++;
	  fprintf(stdout,"Residue is: %s, %d\n", lastSSE.c_str(), sys.getPosition(p).getResidueNumber());
	}
	previousSSE.str("");
      }
      
      lastSSE = sse;
      previousSSE << sse;
      allSSE << sse;
  }
  fprintf(stdout, "FULL SSE: %s\n",allSSE.str().c_str());
  map<string,int>::iterator it;
  for (it = sse_map.begin(); it != sse_map.end();it++){
    fprintf(stdout, "SSE %s %d\n",it->first.c_str(),it->second);
  }
}



bool checkBetaSheet(System &_sys, int _pos){

  Position &pos = _sys.getPosition(_pos);
  if (!(pos.atomExists("N") && pos.atomExists("O"))) return false;

  Atom &posN = pos.getAtom("N");
  Atom &posO = pos.getAtom("O");

  for (uint i = 0; i< _sys.positionSize();i++){
    if (abs(((int)i)-_pos) <= 2) continue;

    Position &pos2 = _sys.getPosition(i);
    if (!(pos2.atomExists("N") && pos2.atomExists("O"))) continue;
    Atom &pos2N = pos2.getAtom("N");
    Atom &pos2O = pos2.getAtom("O");

    int hbonds = 0;
    // Anti-parallel
    // pos == pos2
    double dist = 3.5;
    if (posN.distance(pos2O) < dist){
      hbonds++;
    }
    if (posO.distance(pos2N) < dist){
      hbonds++;
    }

    if (hbonds == 2)  { return true;}
    
    // Anti-parallel
    // Check pos-1 to i+1 AND pos+1 to i-1
    hbonds = 0;
    if (_pos > 0 && _sys.getPosition(_pos-1).atomExists("N") && _sys.getPosition(_pos-1).atomExists("O")){

      if (i < _sys.positionSize()-1 && _sys.getPosition(i+1).atomExists("N") && _sys.getPosition(i+1).atomExists("O")){
	if (_sys.getPosition(_pos-1).getAtom("N").distance(_sys.getPosition(i+1).getAtom("O")) < dist ||
	    _sys.getPosition(_pos-1).getAtom("O").distance(_sys.getPosition(i+1).getAtom("N")) < dist){
	  hbonds++;
	}
      }
    }

    if (_pos < _sys.positionSize()-1 && _sys.getPosition(_pos+1).atomExists("N") && _sys.getPosition(_pos+1).atomExists("O")){

      if (i > 0 && _sys.getPosition(i-1).atomExists("N") && _sys.getPosition(i-1).atomExists("O")){
	if (_sys.getPosition(_pos+1).getAtom("N").distance(_sys.getPosition(i-1).getAtom("O")) < dist ||
	    _sys.getPosition(_pos+1).getAtom("O").distance(_sys.getPosition(i-1).getAtom("N")) < dist){
	  hbonds++;
	}
      }
    }
    if (hbonds == 2) {return true;}


    // Parallel conditions
    // check pos to i-1 AND pos to i+1
    hbonds= 0;
    if (i > 0 && _sys.getPosition(i-1).atomExists("N") && _sys.getPosition(i-1).atomExists("O")){
      
      if (posN.distance(_sys.getPosition(i-1).getAtom("O")) < dist || posO.distance(_sys.getPosition(i-1).getAtom("N")) < dist){
	hbonds++;
      }
    }

    if (i < _sys.positionSize()-1 && _sys.getPosition(i+1).atomExists("N") && _sys.getPosition(i+1).atomExists("O")){
      
      if (posN.distance(_sys.getPosition(i+1).getAtom("O")) < dist || posO.distance(_sys.getPosition(i+1).getAtom("N")) < dist){
	hbonds++;
      }
    }
   
    if (hbonds ==2 ) {return true;}

    // check pos-1 to i AND pos+1 to i
    if (_pos > 0 && _sys.getPosition(_pos-1).atomExists("N") && _sys.getPosition(_pos-1).atomExists("O")){
      
      if (_sys.getPosition(_pos-1).getAtom("N").distance(_sys.getPosition(i).getAtom("O")) < dist ||
	  _sys.getPosition(_pos-1).getAtom("O").distance(_sys.getPosition(i).getAtom("N")) < dist){
	hbonds++;
      }
    }

    if (_pos < _sys.positionSize()-1 && _sys.getPosition(_pos+1).atomExists("N") && _sys.getPosition(_pos+1).atomExists("O")){
      
      if (_sys.getPosition(_pos+1).getAtom("N").distance(_sys.getPosition(i).getAtom("O")) < dist || 
	  _sys.getPosition(_pos+1).getAtom("O").distance(_sys.getPosition(i).getAtom("N")) < dist){
	hbonds++;
      }
    }

    if (hbonds == 2) {return true;}
  }

  return false;
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
  MSLOUT.stream() << "Options:\n"<<OP<<endl;
  return opt;
}
