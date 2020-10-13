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

#include "betaDecider.h"

using namespace std;
using namespace MSL;


// MslOut 
static MslOut MSLOUT("betaDecider");
static SysEnv SYSENV;


int main(int argc, char *argv[]) {

  // Parse commandline options
  Options opt = setupOptions(argc,argv);

  // Read PDB structure
  System sys;
  sys.readPdb(opt.pdb);

  int strandCount = 0;
  vector< pair<int,int> > strands;
  for (uint i = 0; i < sys.positionSize();i++){

    if (checkBetaSheet(sys,i)){
      cout << " SHEET: "<<sys.getPosition(i).getPositionId()<<endl;
      strandCount++;
    } else {

      // Add ending strand to list of strands if more than 4 residues in a row were beta
      if (strandCount > 4){
	strands.push_back(pair<int,int>(i-strandCount,i));
      }

      strandCount = 0;
    }

  }

  
  PyMolVisualization pymol;

  cout << "Number of strands detected: "<<strands.size()<<endl;

  // Print out each strand
  for (uint i = 0; i < strands.size();i++){
    string strandName = MslTools::stringf("Strand-%03d",i+1);


    cout << strandName <<" is between "<<sys.getPosition(strands[i].first).getPositionId()<< " and "<<sys.getPosition(strands[i].second).getPositionId()<<endl;

    // Save for pymol
    stringstream strandSel;
    strandSel << "chain "<<sys.getPosition(strands[i].first).getChainId()<< " and resi "<<sys.getPosition(strands[i].first).getResidueNumber()<<"-"<<sys.getPosition(strands[i].second).getResidueNumber();
    string strandSelStr = strandSel.str();
    string strandColor = "blue";
    pymol.createSelection(strandName,strandSelStr, strandColor);
  }

  ofstream fout;
  fout.open(MslTools::stringf("%s-STRANDS.py", MslTools::getFileName(opt.pdb).c_str()));
  fout << pymol;
  fout.close();

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
    cout << "betaDecider --pdb PDB\n";

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

  return opt;
}



// Check Beta Sheet function
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
    if (posN.distance(pos2O) < 3.25){
      hbonds++;
    }
    if (posO.distance(pos2N) < 3.25){
      hbonds++;
    }

    if (hbonds == 2)  { return true;}
    
    // Anti-parallel
    // Check pos-1 to i+1 AND pos+1 to i-1
    hbonds = 0;
    if (_pos > 0 && _sys.getPosition(_pos-1).atomExists("N") && _sys.getPosition(_pos-1).atomExists("O")){

      if (i < _sys.positionSize()-1 && _sys.getPosition(i+1).atomExists("N") && _sys.getPosition(i+1).atomExists("O")){
	if (_sys.getPosition(_pos-1).getAtom("N").distance(_sys.getPosition(i+1).getAtom("O")) < 3.25 ||
	    _sys.getPosition(_pos-1).getAtom("O").distance(_sys.getPosition(i+1).getAtom("N")) < 3.25){
	  hbonds++;
	}
      }
    }

    if (_pos < _sys.positionSize()-1 && _sys.getPosition(_pos+1).atomExists("N") && _sys.getPosition(_pos+1).atomExists("O")){

      if (i > 0 && _sys.getPosition(i-1).atomExists("N") && _sys.getPosition(i-1).atomExists("O")){
	if (_sys.getPosition(_pos+1).getAtom("N").distance(_sys.getPosition(i-1).getAtom("O")) < 3.25 ||
	    _sys.getPosition(_pos+1).getAtom("O").distance(_sys.getPosition(i-1).getAtom("N")) < 3.25){
	  hbonds++;
	}
      }
    }
    if (hbonds == 2) {return true;}


    // Parallel conditions
    // check pos to i-1 AND pos to i+1
    hbonds= 0;
    if (i > 0 && _sys.getPosition(i-1).atomExists("N") && _sys.getPosition(i-1).atomExists("O")){
      
      if (posN.distance(_sys.getPosition(i-1).getAtom("O")) < 3.25 || posO.distance(_sys.getPosition(i-1).getAtom("N")) < 3.25){
	hbonds++;
      }
    }

    if (i < _sys.positionSize()-1 && _sys.getPosition(i+1).atomExists("N") && _sys.getPosition(i+1).atomExists("O")){
      
      if (posN.distance(_sys.getPosition(i+1).getAtom("O")) < 3.25 || posO.distance(_sys.getPosition(i+1).getAtom("N")) < 3.25){
	hbonds++;
      }
    }
   
    if (hbonds ==2 ) {return true;}

    // check pos-1 to i AND pos+1 to i
    if (_pos > 0 && _sys.getPosition(_pos-1).atomExists("N") && _sys.getPosition(_pos-1).atomExists("O")){
      
      if (_sys.getPosition(_pos-1).getAtom("N").distance(_sys.getPosition(i).getAtom("O")) < 3.25 ||
	  _sys.getPosition(_pos-1).getAtom("O").distance(_sys.getPosition(i).getAtom("N")) < 3.25){
	hbonds++;
      }
    }

    if (_pos < _sys.positionSize()-1 && _sys.getPosition(_pos+1).atomExists("N") && _sys.getPosition(_pos+1).atomExists("O")){
      
      if (_sys.getPosition(_pos+1).getAtom("N").distance(_sys.getPosition(i).getAtom("O")) < 3.25 || 
	  _sys.getPosition(_pos+1).getAtom("O").distance(_sys.getPosition(i).getAtom("N")) < 3.25){
	hbonds++;
      }
    }

    if (hbonds == 2) {return true;}
  }

  return false;
}
