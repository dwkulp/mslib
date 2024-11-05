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
#include "Frame.h"
#include "AtomPointerVector.h"
#include "PyMolVisualization.h"

#include "betaAnalysis.h"

using namespace std;
using namespace MSL;


// MslOut for nice output
static MslOut MSLOUT("betaDecider");
static SysEnv SYSENV;


int main(int argc, char *argv[]) {

  // Parse commandline options
  Options opt = setupOptions(argc,argv);

  // Read PDB structure
  System sys;
  sys.readPdb(opt.pdb);

  // Count the number of strand residues
  int strandCount = 0;

  // Store list of the position for start, end of each strand
  vector< pair<int,int> > strands;

  // Loop over all positions
  for (uint i = 0; i < sys.positionSize();i++){

    // Check if this position is in a beta sheet
    if (i > 0 && i < sys.positionSize()-2){
      double phi = PhiPsiStatistics::getPhi(sys.getPosition(i-1).getCurrentIdentity(),sys.getPosition(i).getCurrentIdentity());
      double psi = PhiPsiStatistics::getPsi(sys.getPosition(i).getCurrentIdentity(),sys.getPosition(i+1).getCurrentIdentity());
      if (phi < -35 &&
	  phi > -90 &&
	  psi > -70 &&
	  psi < 0){
	continue;
      }
    }

    if (checkBetaSheet(sys,i)){

      // Print out only in debug mode
      MSLOUT.debug() << " SHEET: "<<sys.getPosition(i).getPositionId()<<endl;


      // Increment Strand Count
      strandCount++;


    } else {

      // This residue is not a strand, so could be end of strand

      // Add ending strand to list of strands if more than 4 residues in a row were beta
      if (strandCount > 4){
	strands.push_back(pair<int,int>(i-strandCount,i));
      }

      // Reset strand count
      strandCount = 0;
    }

  }



  // Output strand defintions to screen and make a PyMOL script containing selections for each strand
  PyMolVisualization pymol;

  cout << "Number of strands detected: "<<strands.size()<<endl;

  // Print out each strand
  for (uint i = 0; i < strands.size();i++){

    // Give strand a name
    string strandName = MslTools::stringf("Strand-%03d",i+1);

    // Print out strand defition
    cout << strandName <<" is between "<<sys.getPosition(strands[i].first).getPositionId()<< " and "<<sys.getPosition(strands[i].second).getPositionId()<<endl;

    // Save for pymol
    stringstream strandSel;
    strandSel << "chain "<<sys.getPosition(strands[i].first).getChainId()<< " and resi "<<sys.getPosition(strands[i].first).getResidueNumber()<<"-"<<sys.getPosition(strands[i].second).getResidueNumber();
    string strandSelStr = strandSel.str();
    string strandColor = "blue";
    pymol.createSelection(strandName,strandSelStr, strandColor);
  }





  // Create an adjacency list
  vector<vector<int> > adList;
  for (uint i = 0; i < strands.size();i++){
    
    vector<int> neighbors_of_strand_i;
    for (uint j = 0; j < strands.size();j++){
      if (i == j) continue;
      if (neighboringStrand(sys,strands[i],strands[j])){
	neighbors_of_strand_i.push_back(j);
      }
    }
    adList.push_back(neighbors_of_strand_i);
  }

  for (uint i = 0; i < adList.size();i++){
    cout << "Strand "<<MslTools::stringf("%-3d",i+1)<<" adjacency list: ";
    for (uint j = 0; j < adList[i].size();j++){
      cout << MslTools::stringf("%-3d",adList[i][j]+1) << " ";
    }
    cout << endl;
  }



  // Get Frames for each set of 4 CA atoms
  vector<Frame> frames;
  vector<int> strand_index;
  AtomPointerVector allCAs;
  for (uint i  = 0; i < strands.size();i++){
    for (uint j = strands[i].first; j<=strands[i].second-4;j+=2){
      AtomPointerVector apv;
      apv.push_back(&sys.getPosition(j).getAtom("CA"));
      apv.push_back(&sys.getPosition(j+1).getAtom("CA"));
      apv.push_back(&sys.getPosition(j+2).getAtom("CA"));
      apv.push_back(&sys.getPosition(j+3).getAtom("CA"));

      Frame fr;
      fr.setName(MslTools::stringf("Fr%02d", frames.size()));
      fr.computeFrameFromPCA(apv);
      frames.push_back(fr);
      strand_index.push_back(i);
    }
    for (uint j = strands[i].first; j <= strands[i].first;j++){
      allCAs.push_back(&sys.getPosition(j).getAtom("CA"));
    }
  }

  Frame bFrame;
  bFrame.setName("bFrame");
  bFrame.computeFrameFromPCA(allCAs);

  // Write out PyMOL Script in PyMol use "run PDB-STRANDS.py"
  ofstream fout;
  fout.open(MslTools::stringf("%s-STRANDS.py", MslTools::getFileName(opt.pdb).c_str()));
  fout << pymol;
  fout << bFrame.toString();

  for (uint f = 0; f<frames.size();f++){
    double dist = bFrame.distanceToFrame(frames[f]);
    double ang  = bFrame.anglesBetweenFrame(frames[f])[0][0];
    
    

    cout << MslTools::stringf("%d, %d, %8.2f, %8.2f\n",strand_index[f],f,dist,ang);
    fout << frames[f].toString();
  }  







  fout.close();
  
}



// Function to help setup input options
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
    if (posN.distance(pos2O) < 3.5){
      hbonds++;
    }
    if (posO.distance(pos2N) < 3.5){
      hbonds++;
    }

    if (hbonds == 2)  { return true;}
    
    // Anti-parallel
    // Check pos-1 to i+1 AND pos+1 to i-1
    hbonds = 0;
    if (_pos > 0 && _sys.getPosition(_pos-1).atomExists("N") && _sys.getPosition(_pos-1).atomExists("O")){

      if (i < _sys.positionSize()-1 && _sys.getPosition(i+1).atomExists("N") && _sys.getPosition(i+1).atomExists("O")){
	if (_sys.getPosition(_pos-1).getAtom("N").distance(_sys.getPosition(i+1).getAtom("O")) < 3.5 ||
	    _sys.getPosition(_pos-1).getAtom("O").distance(_sys.getPosition(i+1).getAtom("N")) < 3.5){
	  hbonds++;
	}
      }
    }

    if (_pos < _sys.positionSize()-1 && _sys.getPosition(_pos+1).atomExists("N") && _sys.getPosition(_pos+1).atomExists("O")){

      if (i > 0 && _sys.getPosition(i-1).atomExists("N") && _sys.getPosition(i-1).atomExists("O")){
	if (_sys.getPosition(_pos+1).getAtom("N").distance(_sys.getPosition(i-1).getAtom("O")) < 3.5 ||
	    _sys.getPosition(_pos+1).getAtom("O").distance(_sys.getPosition(i-1).getAtom("N")) < 3.5){
	  hbonds++;
	}
      }
    }
    if (hbonds == 2) {return true;}


    // Parallel conditions
    // check pos to i-1 AND pos to i+1
    hbonds= 0;
    if (i > 0 && _sys.getPosition(i-1).atomExists("N") && _sys.getPosition(i-1).atomExists("O")){
      
      if (posN.distance(_sys.getPosition(i-1).getAtom("O")) < 3.5 || posO.distance(_sys.getPosition(i-1).getAtom("N")) < 3.5){
	hbonds++;
      }
    }

    if (i < _sys.positionSize()-1 && _sys.getPosition(i+1).atomExists("N") && _sys.getPosition(i+1).atomExists("O")){
      
      if (posN.distance(_sys.getPosition(i+1).getAtom("O")) < 3.5 || posO.distance(_sys.getPosition(i+1).getAtom("N")) < 3.5){
	hbonds++;
      }
    }
   
    if (hbonds ==2 ) {return true;}

    // check pos-1 to i AND pos+1 to i
    if (_pos > 0 && _sys.getPosition(_pos-1).atomExists("N") && _sys.getPosition(_pos-1).atomExists("O")){
      
      if (_sys.getPosition(_pos-1).getAtom("N").distance(_sys.getPosition(i).getAtom("O")) < 3.5 ||
	  _sys.getPosition(_pos-1).getAtom("O").distance(_sys.getPosition(i).getAtom("N")) < 3.5){
	hbonds++;
      }
    }

    if (_pos < _sys.positionSize()-1 && _sys.getPosition(_pos+1).atomExists("N") && _sys.getPosition(_pos+1).atomExists("O")){
      
      if (_sys.getPosition(_pos+1).getAtom("N").distance(_sys.getPosition(i).getAtom("O")) < 3.5 || 
	  _sys.getPosition(_pos+1).getAtom("O").distance(_sys.getPosition(i).getAtom("N")) < 3.5){
	hbonds++;
      }
    }

    if (hbonds == 2) {return true;}
  }

  return false;
}


bool neighboringStrand(System &_sys, pair<int,int> &_s1, pair<int,int> &_s2){

  
  int hbonds=0;
  for (uint i = _s1.first; i <= _s1.second;i++){
    for (uint j = _s2.first; j <= _s2.second;j++){
      if (_sys.getPosition(i).atomExists("N") && _sys.getPosition(j).atomExists("O")){
	if (_sys.getPosition(i).getAtom("N").distance(_sys.getPosition(j).getAtom("O")) < 3.5){
	  hbonds++;
	}
      }
      if (_sys.getPosition(i).atomExists("O") && _sys.getPosition(j).atomExists("N")){
	if (_sys.getPosition(i).getAtom("O").distance(_sys.getPosition(j).getAtom("N")) < 3.5){
	  hbonds++;
	}
      }
    }
  }

  if (hbonds > 1){
    return true;
  }

  return false;

}
