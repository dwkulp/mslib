
#include "System.h"
#include "OptionParser.h"
#include "SasaCalculator.h"
#include "PhiPsiStatistics.h"
#include "PyMolVisualization.h"
#include "MslTools.h"
#include "MslOut.h"
#include "AtomSelection.h"
#include "getExposedPositions.h"
#include <iostream>
#include <fstream>

using namespace MSL;
using namespace std;

static MslOut MSLOUT("getExposedPositions");

int main(int argc, char *argv[]) {
    Options opt = setupOptions(argc, argv);

    System sys;
    if (!sys.readPdb(opt.pdbFile)) {
        cerr << "Unable to read PDB file: " << opt.pdbFile << endl;
        return 1;
    }

    vector<subSeq> exposedPositions = getExposedPositions(sys, opt);
    ofstream outFile(MslTools::getFileName(opt.pdbFile) + "_exposed_positions.txt");
    for (size_t i = 0; i < exposedPositions[0].positions.size(); ++i) {
        Position &pos = sys.getPosition(exposedPositions[0].positions[i]);
        outFile << pos.getChainId() << pos.getResidueNumber() << ",";
    }
    outFile << endl;
    outFile.close();

    // Write exposed positions to JSON file
    ofstream jsonFile(MslTools::getFileName(opt.pdbFile) + "_exposed_positions.json");
    jsonFile << "{\"" << MslTools::getFileName(opt.pdbFile) << "\": { \""<<sys.getPosition(exposedPositions[0].positions[0]).getChainId()<<"\": [";
    for (size_t i = 0; i < exposedPositions[0].positions.size(); ++i) {
      Position &pos = sys.getPosition(exposedPositions[0].positions[i]);
      jsonFile << pos.getResidueNumber();
      if (i < exposedPositions[0].positions.size() - 1) {
        jsonFile << ",";
      }
    }
    jsonFile << "]}}" << endl;
    jsonFile.close();

    ofstream outFileBuried(MslTools::getFileName(opt.pdbFile) + "_buried_positions.txt");
    for (size_t i = 0; i < exposedPositions[1].positions.size(); ++i) {
        Position &pos = sys.getPosition(exposedPositions[1].positions[i]);
        outFileBuried << pos.getChainId() << pos.getResidueNumber() << ",";
    }
    outFileBuried << endl;
    outFileBuried.close();
    // Write buried positions to JSON file
    ofstream jsonFileBuried(MslTools::getFileName(opt.pdbFile) + "_buried_positions.json");
    jsonFileBuried << "{\"" << MslTools::getFileName(opt.pdbFile) << "\": { \""<<sys.getPosition(exposedPositions[1].positions[0]).getChainId()<<"\": [";
    for (size_t i = 0; i < exposedPositions[1].positions.size(); ++i) {
      Position &pos = sys.getPosition(exposedPositions[1].positions[i]);
      jsonFileBuried << pos.getResidueNumber();
      if (i < exposedPositions[1].positions.size() - 1) {
      jsonFileBuried << ",";
      }
    }
    jsonFileBuried << "]}}" << endl;
    jsonFileBuried.close();
    return 0;
}

Options setupOptions(int argc, char *argv[]) {
    Options opt;
    OptionParser OP;

    OP.readArgv(argc, argv);
    OP.setRequired(opt.required);
    OP.setRequired(opt.optional);
    

    opt.pdbFile = OP.getString("pdb");
    opt.percentSasa = OP.getDouble("percentSasa");
    opt.selectPositions = OP.getString("selectPositions");
    opt.pymol = OP.getBool("pymol");

    return opt;
}

vector<subSeq> getExposedPositions(System &_sys, Options &_opt){

  PyMolVisualization pymol;

  vector<subSeq> results(2);

     /*
	SASA reference:
	Protein Engineering vol.15 no.8 pp.659–667, 2002
	Quantifying the accessible surface area of protein residues in their local environment
	Uttamkumar Samanta Ranjit P.Bahadur and  Pinak Chakrabarti
      */
  map<string, double> refSasa;
  refSasa["G"] = 83.91;
  refSasa["A"] = 116.40;
  refSasa["S"] = 125.68;
  refSasa["C"] = 141.48;
  refSasa["P"] = 144.80;
  refSasa["T"] = 148.06;
  refSasa["D"] = 155.37;
  refSasa["V"] = 162.24;
  refSasa["N"] = 168.87;
  refSasa["E"] = 187.16;
  refSasa["Q"] = 189.17;
  refSasa["I"] = 189.95;
  refSasa["L"] = 197.99;
  refSasa["H"] = 198.51;
  refSasa["K"] = 207.49;
  refSasa["M"] = 210.55;
  refSasa["F"] = 223.29;
  refSasa["Y"] = 238.30;
  refSasa["R"] = 249.26;
  refSasa["W"] = 265.42;

  // Identify Exposed positions
  SasaCalculator scalc(_sys.getAtomPointers());
  scalc.calcSasa();

  string selpos = "all";
  if (_opt.selectPositions != ""){
    selpos = _opt.selectPositions;
  }
  AtomSelection atsel(_sys.getAtomPointers());
  AtomPointerVector &av = atsel.select("resurface,"+selpos);
  MSLOUT.stream() << "Number of selected atoms: "<<av.size()<<" "<<selpos<<endl;
  
  stringstream selstr;
  selstr << "resi ";
  for (uint p = 1; p < _sys.positionSize()-1;p++){
    if (!_sys.getPosition(p).atomExists("CA")) continue;
    if (!_sys.getPosition(p).getAtom("CA").getSelectionFlag("resurface")) continue;

    double normSasa = scalc.getResidueSasa(_sys.getPosition(p).getPositionId()) / refSasa[MslTools::getOneLetterCode(_sys.getPosition(p).getResidueName())];
    if (normSasa > 1.0){
      normSasa = 1.0;
    }
    //double phi = PhiPsiStatistics::getPhi(_sys.getPosition(p-1).getCurrentIdentity(),_sys.getPosition(p).getCurrentIdentity());
    //double psi = PhiPsiStatistics::getPsi(_sys.getPosition(p).getCurrentIdentity(),_sys.getPosition(p+1).getCurrentIdentity());


    // Exposed if > 40% ?
    //string msg   = MslTools::stringf("Position %7s %8.2f %8.2f %8.2f %8.2f %s",_sys.getPosition(p).getPositionId().c_str(),scalc.getResidueSasa(_sys.getPosition(p).getPositionId()),normSasa,phi,psi,sse.c_str());      
    if (normSasa > _opt.percentSasa) {
    //msg += MslTools::stringf(" **** ");
      
      results[0].seq += MslTools::getOneLetterCode(_sys.getPosition(p).getResidueName());
      results[0].positions.push_back(p);
      //MSLOUT.stream() << "Exposed position "<<_sys.getPosition(p).getPositionId()<<" "<<MslTools::getOneLetterCode(_sys.getPosition(p).getResidueName())<< " "<<result.sse<<endl;
      selstr << _sys.getPosition(p).getResidueNumber()<<"+";


    } else {
        results[1].seq += MslTools::getOneLetterCode(_sys.getPosition(p).getResidueName());
        results[1].positions.push_back(p);
    }
    //MSLOUT.stream() << MslTools::stringf("%s\n",msg.c_str());

  } // END FOR POSITIONS
  if (_opt.pymol){
    string selname = "exposed";
    string sel = selstr.str();
      pymol.createSelection(selname,sel);
    cout  << pymol.toString()<<endl;
    }       
  
  return results;
}