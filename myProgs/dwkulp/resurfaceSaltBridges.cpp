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

#include "resurfaceSaltBridges.h"

using namespace std;
using namespace MSL;


// MslOut 
static MslOut MSLOUT("resurfaceSaltBridges");
static SysEnv SYSENV;

struct subSeq {
  string seq;
  vector<int> positions;
  string sse;
  subSeq(){
    seq = "";
    sse = "";
  }
  subSeq(const subSeq &_rhs){
    seq       = _rhs.seq;
    positions = _rhs.positions;
    sse       = _rhs.sse;
  }
  void operator=(const subSeq &_rhs){
    seq       = _rhs.seq;
    positions = _rhs.positions;
    sse       = _rhs.sse;
  }
};

struct saltBridgeResult {
  double score;
  vector<pair<string,string> > mutations;
  string sbId;
  
  saltBridgeResult(){
    score = 1.0;
    sbId  = "";
  }

  saltBridgeResult(const saltBridgeResult &_rhs){
    score       = _rhs.score;
    mutations   = _rhs.mutations;
    sbId        = _rhs.sbId;
  }
  void operator=(const saltBridgeResult &_rhs){
    score       = _rhs.score;
    mutations   = _rhs.mutations;
    sbId        = _rhs.sbId;
  }

};

struct compareScores {
  bool operator() (pair<saltBridgeResult, subSeq> &lhs,pair<saltBridgeResult, subSeq> &rhs){
    return (lhs.first.score < rhs.first.score);
  }
};

// Local utility functions
subSeq getExposedPositions(System &_sys, Options &_opt);
bool checkBetaSheet(System &_sys, int _pos);
void readSaltBridgeData(string _sbFile, map<string,map<string,map<int, map<string, map<string,double> > > > > &_sb_table);
saltBridgeResult scoreSaltBridgePositions(subSeq &_seq, map<string,map<string,map<int, map<string, map<string,double> > > > > &_sb_table, System &_sys, bool verbose=false);
void runSaltBridgeMC(subSeq &_exposedPositions,
		     map<string,map<string,map<int, map<string,map<string, double> > > > > &_sbScoreTable,
		     System &_sys,
		     Options &_opt,
		     std::priority_queue< std::pair<saltBridgeResult,subSeq>, std::vector< std::pair<saltBridgeResult,subSeq> >, compareScores> &_mcData);

void runModelSaltBridges(System &_sys,
			 Options &_opt,
			 std::priority_queue< std::pair<saltBridgeResult,subSeq>, std::vector< std::pair<saltBridgeResult,subSeq> >, compareScores> &_mcData);
string mutate(System &_sys, string _mutationId, Options &_opt);
pair<string,string> getPosIdAndMutation(System &_sys, string _mutationId, Options &_opt);

int main(int argc, char *argv[]) {

  // Parse commandline options
  Options opt = setupOptions(argc,argv);

  // Read PDB structure
  System sys;
  sys.readPdb(opt.pdb);

  // Discover exposed positions
  subSeq exposedPositions = getExposedPositions(sys,opt);
  
  // Read Salt Bridge Propensity Table
  // ['BasicAA']['basicSSE'][acidicPos]['AcidicAA']['AcidicSS'] = -ln(ratio)
  map<string,map<string,map<int, map<string,map<string, double> > > > > sbScoreTable;
  readSaltBridgeData(opt.sb_prop_table,sbScoreTable);

  // Score exposed positions on protein
  saltBridgeResult wt_sbresult;
  if (opt.scoreOnly){
    wt_sbresult = scoreSaltBridgePositions(exposedPositions,sbScoreTable,sys,true);
    fprintf(stdout, "SCORE: %8.3f\n",wt_sbresult.score);
    for (uint i = 0; i < wt_sbresult.mutations.size();i++){
      fprintf(stdout,"\t%s %s\n",wt_sbresult.mutations[i].first.c_str(),wt_sbresult.mutations[i].second.c_str());
    }
    exit(0);
  } else {
    wt_sbresult = scoreSaltBridgePositions(exposedPositions,sbScoreTable,sys);
    MSLOUT.stream() << "WT SCORE: "<<wt_sbresult.score<<endl;
  }

  // Priority queue to store top scoring solutions
  std::priority_queue< std::pair<saltBridgeResult,subSeq>, std::vector< std::pair<saltBridgeResult,subSeq> >, compareScores> mcData;

  // Run MC for SaltBridges
  runSaltBridgeMC(exposedPositions, sbScoreTable, sys, opt, mcData);


  // Model mutations
  runModelSaltBridges(sys,opt,mcData);


  MSLOUT.stream() << "Done."<<endl;
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
  opt.sb_prop_table = OP.getString("sb_prop_table");
  if (OP.fail()){
    opt.sb_prop_table = SYSENV.getEnv("MSL_SB_PROP_TABLE"); 
    cout << "SB_PROP_TABLE: "<<opt.sb_prop_table<<endl;
    if (opt.sb_prop_table == "UNDEF"){
      std::cerr << "ERROR 1111 no sb_prop_table specified."<<std::endl;	
      exit(1111);
    } else {
      if (MslTools::fileExists(opt.sb_prop_table)){
	cerr << "WARNING sb_prop_table defaulted to: "<<opt.sb_prop_table<<endl;
      } else {
	std::cerr << "ERROR 1111 no sb_prop_table specified and couldn't find default one: "<<opt.sb_prop_table<<std::endl;	
	exit(1111);	
      }
      
    }

  }

  opt.topfile = OP.getString("topfile");
  if (OP.fail()){
    opt.topfile = SYSENV.getEnv("MSL_CHARMM_TOP");
    cerr << "WARNING: charmmtopfile not specified, using " << opt.topfile << "\n";
  }
  if (!MslTools::fileExists(opt.topfile)){
    cerr << "ERROR 1111 CHARMM TOPFILE DOESN'T EXIST: "<<opt.topfile<<"."<<endl;
    exit(1111);
  }

  opt.parfile = OP.getString("parfile");
  if (OP.fail()){
    opt.parfile = SYSENV.getEnv("MSL_CHARMM_PAR");
    cerr <<"WARNING charmmparfile not specified, using " << opt.parfile << "\n";

  }
  if (!MslTools::fileExists(opt.parfile)){
    cerr << "ERROR 1111 CHARMM parFILE DOESN'T EXIST: "<<opt.parfile<<endl;
    exit(1111);
  }
  opt.rotlib = OP.getString("rotlib");
  if (OP.fail()){
    opt.rotlib = SYSENV.getEnv("MSL_ROTLIB");
    cerr << "WARNING rotlib not specified, using " << opt.rotlib << "\n";
  }

  if (!MslTools::fileExists(opt.rotlib)){
    cerr << "ERROR 1111 rotamer library DOESN'T EXIST: "<<opt.rotlib<<endl;
    exit(1111);
  }

  opt.numStructuralModels = OP.getInt("numStructuralModels");
  if (OP.fail()){
    opt.numStructuralModels = 10;
  }
  opt.numSequenceModels   = OP.getInt("numSequenceModels");
  if (OP.fail()){
    opt.numSequenceModels = 30;
  }
  opt.numMCcycles  = OP.getInt("numMCcycles");
  if (OP.fail()){
    opt.numMCcycles = 100;
  }
  opt.startMCtemp  = OP.getDouble("startMCtemp");
  if (OP.fail()){
    opt.startMCtemp = 100.0;
  }
  opt.endMCtemp    = OP.getDouble("endMCtemp");
  if (OP.fail()){
    opt.endMCtemp = 1.0;
  }

  opt.scoreOnly = OP.getBool("scoreOnly");

  opt.selectPositions  = OP.getString("selectPositions");
  if (OP.fail()){
    opt.selectPositions = "";
  }

  opt.percentSasa  = OP.getDouble("percentSasa");
  if (OP.fail()){
    opt.percentSasa = 0.4;
  }
  
  opt.fasta_only = OP.getBool("fasta");
  opt.appendName = OP.getString("appendName");
  if (OP.fail()){
    opt.appendName = "";
  } else {
    // Remove \ if first character (needed for strings with "--MSG--" for instance)
    if (opt.appendName.substr(0,1) == "\\"){
      opt.appendName = opt.appendName.substr(1, opt.appendName.size() - 1);
    }
  }
  MSLOUT.stream() << "Options:\n"<<OP<<endl;
  return opt;
}


void readSaltBridgeData(string _sbFile, map<string,map<string,map<int, map<string, map<string,double> > > > > &_sb_table){


  ifstream in;
  in.open(_sbFile.c_str());
  if(!in.is_open()) {
    cerr << "unable to open " << _sbFile << endl;
    exit(0);
  }
  string line;
  while(getline(in,line)) {
    // Skip comment lines
    if (line[0] == '#') continue;

    // Skip blank lines
    if (MslTools::trim(line).size() == 0) continue;


    vector<string> tokens = MslTools::tokenize(line);
    if (tokens.size() != 8) {
      cerr << "ERROR 3444 in salt bridge table, expected 8 tokens got: "<<tokens.size()<<endl;
      cerr << "\tLINE:"<<line<<endl;
      exit(3444);
    }

    _sb_table[tokens[0]][tokens[2]][MslTools::toInt(tokens[4])][tokens[1]][tokens[3]] = -log(MslTools::toDouble(tokens[7]));
    cout << "SB TABLE: "<<tokens[0]<<","<<tokens[2]<<","<<MslTools::toInt(tokens[4])<<","<<tokens[1]<<","<<tokens[3]<<","<< -log(MslTools::toDouble(tokens[7]))<<endl;
  }
  in.close();

}


subSeq getExposedPositions(System &_sys, Options &_opt){

  PyMolVisualization pymol;

  subSeq result;

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
    //if (!_sys.getPosition(p).atomExists("CA")) continue;
    if (!_sys.getPosition(p).getAtom("CA").getSelectionFlag("resurface")) continue;

    double normSasa = scalc.getResidueSasa(_sys.getPosition(p).getPositionId()) / refSasa[MslTools::getOneLetterCode(_sys.getPosition(p).getResidueName())];
    if (normSasa > 1.0){
      normSasa = 1.0;
    }
    double phi = PhiPsiStatistics::getPhi(_sys.getPosition(p-1).getCurrentIdentity(),_sys.getPosition(p).getCurrentIdentity());
    double psi = PhiPsiStatistics::getPsi(_sys.getPosition(p).getCurrentIdentity(),_sys.getPosition(p+1).getCurrentIdentity());


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
	bool isSheet = checkBetaSheet(_sys, p);
	if (isSheet){
	  sse = "E";
	} else {
	  sse = "O";
	}
      }


    // Exposed if > 40% ?
    string msg   = MslTools::stringf("Position %7s %8.2f %8.2f %8.2f %8.2f %s",_sys.getPosition(p).getPositionId().c_str(),scalc.getResidueSasa(_sys.getPosition(p).getPositionId()),normSasa,phi,psi,sse.c_str());      
    if (normSasa > _opt.percentSasa) {
       msg += MslTools::stringf(" **** ");
      result.sse += sse;
      result.seq += MslTools::getOneLetterCode(_sys.getPosition(p).getResidueName());
      result.positions.push_back(p);
      //MSLOUT.stream() << "Exposed position "<<_sys.getPosition(p).getPositionId()<<" "<<MslTools::getOneLetterCode(_sys.getPosition(p).getResidueName())<< " "<<result.sse<<endl;
      selstr << _sys.getPosition(p).getResidueNumber()<<"+";


    } // IF EXPOSED POSITION
    MSLOUT.stream() << MslTools::stringf("%s\n",msg.c_str());

  } // END FOR POSITIONS
  string selname = "exposed";
  string sel = selstr.str();
  pymol.createSelection(selname,sel);
  cout << pymol.toString()<<endl;
  return result;
}


saltBridgeResult scoreSaltBridgePositions(subSeq &_seqDat, map<string,map<string,map<int, map<string,map<string, double> > > > > &_sb_table, System &_sys, bool verbose) {

  if (verbose){
    MSLOUT.stream() << "SEQ: "<<_seqDat.seq<<endl;
    MSLOUT.stream() << "SSE: "<<_seqDat.sse<<endl;
  }
  stringstream ss;
  ss << "POS: ";
  for (uint i = 0 ; i < _seqDat.positions.size();i++){
    ss << _seqDat.positions[i] << " ";
  }
  MSLOUT.debug()<< ss.str()<<endl;

  saltBridgeResult sbr;
  sbr.score = 0.0;
  for (uint i = 0; i < _seqDat.seq.size();i++){

    MSLOUT.debug() << "Scoring position "<<i<<" it is "<<_seqDat.seq[i]<<" "<<_seqDat.sse[i]<<" "<<_seqDat.positions[i]<<endl;
    // Only score on the basic residue
    if (_seqDat.seq[i] != 'R' && _seqDat.seq[i] != 'K' && _seqDat.seq[i] != 'H') continue;


    string basic_aa  = MslTools::getThreeLetterCode(_seqDat.seq.substr(i,1));
    string basic_sse = _seqDat.sse.substr(i,1);

    map<string,map<string,map<int, map<string,map<string, double> > > > >::iterator basicAAit;
    basicAAit = _sb_table.find(basic_aa);
    if (basicAAit == _sb_table.end()){
      throw MslNotFoundException(MslTools::stringf("ERROR Basic residue AA: %s not found in salt-bridge table", basic_aa.c_str()));
    }

    map<string,map<int, map<string,map<string, double> > > >::iterator basicSSEit;
    basicSSEit = basicAAit->second.find(basic_sse);
    if (basicSSEit == basicAAit->second.end()){
      MSLOUT.debug() <<MslTools::stringf("Basic residue SSE: %s-%s not found in salt-bridge table", basic_aa.c_str(),basic_sse.c_str())<<endl;
      continue;
    }
    MSLOUT.debug() << "BasicAA and BasicSSE satisified"<<endl;

    // Loop over all local positions..
    for (int pos = -4; pos <= 4; pos++){
      map<int, map<string,map<string, double> > >::iterator posIt;
      posIt = basicSSEit->second.find(pos);
      if (posIt == basicSSEit->second.end()){
	MSLOUT.debug() << "No data for position "<<pos<<" residue "<<basic_aa<<" "<<basic_sse<<endl;
	continue;
      }
      MSLOUT.debug() << "Separation satisified"<<endl;
      int posToCheck = _seqDat.positions[i] + pos;

      // Iterate over positions in _seqDat.positions looking for 'posToCheck' 
      for (int pos2 = -4; pos2 <= 4; pos2++){

	if (i+pos2 < 0) continue;
	if (i+pos2 > _seqDat.positions.size()) continue;

	if (_seqDat.positions[i+pos2] == posToCheck) {
	  map<string,map<string, double> >::iterator acidAAit;
	  acidAAit = posIt->second.find(MslTools::getThreeLetterCode(_seqDat.seq.substr(i+pos2,1)));
	  if (acidAAit == posIt->second.end()){
	    MSLOUT.debug() << "Position "<<i+pos2<<" doesn't have an acidic residue."<<endl;
	    continue;
	  }
	  MSLOUT.debug() << "AcidicAA satisified"<<endl;
	  // LAST ONE!
	  map<string, double>::iterator acidSSEit;
	  acidSSEit = acidAAit->second.find(_seqDat.sse.substr(i+pos2,1));
	  if (acidSSEit == acidAAit->second.end()){
	    MSLOUT.debug() << "Position "<<i+pos2<<" has acidic residue: "<<_seqDat.seq[i+pos2]<<" but doesn't have an acidic residue sse: "<<_seqDat.sse[i+pos2]<<endl;
	    continue;
	  }

	  // WE HAVE IT!
	  if (verbose){
	    MSLOUT.stream() << " HAS A SALT BRIDGE: "<<_sys.getPosition(_seqDat.positions[i]).getCurrentIdentity().getIdentityId()<<" to "<<basic_aa<<" and "<<_sys.getPosition(_seqDat.positions[i+pos2]).getCurrentIdentity().getIdentityId()<<" to "<<acidAAit->first <<" with value: "<<acidSSEit->second<<endl;
	  }
	  sbr.score += acidSSEit->second;
	  sbr.mutations.push_back(pair<string,string>(_sys.getPosition(_seqDat.positions[i]).getCurrentIdentity().getIdentityId() + "," + basic_aa,_sys.getPosition(_seqDat.positions[i+pos2]).getCurrentIdentity().getIdentityId()+","+acidAAit->first)) ;
	  sbr.sbId += _sys.getPosition(_seqDat.positions[i]).getCurrentIdentity().getIdentityId() + ":" + basic_aa+" = "+_sys.getPosition(_seqDat.positions[i+pos2]).getCurrentIdentity().getIdentityId()+":"+acidAAit->first+"\n";
	  
	  break;
	} // IF posToCheck

      } // END FOR POS2
    } // END FOR POS
  } // END _seqDat.seq

  return sbr;
  
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


void runSaltBridgeMC(
		     subSeq &_exposedPositions,
		     map<string,map<string,map<int, map<string,map<string, double> > > > > &_sbScoreTable,
		     System &_sys, 
		     Options &_opt,
		     std::priority_queue< std::pair<saltBridgeResult,subSeq>, std::vector< std::pair<saltBridgeResult,subSeq> >, compareScores> &_mcData){

  // WT score
  saltBridgeResult wt_sbresult = scoreSaltBridgePositions(_exposedPositions,_sbScoreTable,_sys);

  // Map to insure we don't fill priority queue with duplicates
  map<string,bool> mcDataMap;

  string possibleMutations = "KREDH";

  RandomNumberGenerator rng;
  rng.setTimeBasedSeed();
  MSLOUT.stream() << "RNG SEED: "<<rng.getSeed()<<endl;

  //MonteCarloManager MCMngr(_startingTemperature, _endingTemperature,_scheduleCycles, _scheduleShape, _maxRejectionsNumber, _convergedSteps, _convergedE);
  MonteCarloManager MCMngr(_opt.startMCtemp,_opt.endMCtemp,_opt.numMCcycles, MonteCarloManager::EXPONENTIAL,100,0,0.0);
  MCMngr.setRandomNumberGenerator(&rng);
  MCMngr.setEner(wt_sbresult.score);
  
  subSeq current   = _exposedPositions;
  string wtSeq     = _exposedPositions.seq;

  // TODO add options to Options opt; to take care of MC-related options..... get rid of hard-coded values.
  saltBridgeResult sbResult = scoreSaltBridgePositions(current,_sbScoreTable,_sys);
  mcDataMap[sbResult.sbId] = true;
  _mcData.push(pair<saltBridgeResult,subSeq>(sbResult,current));

  int counter = 0;
  int total = _opt.numMCcycles*10;
  while (!MCMngr.getComplete()) {
    if (counter > total){
      break;
    }
    
    counter++;
    
    // Make change
    int mutIndex = rng.getRandomInt(0,possibleMutations.size()-1);
    int seqPosIndex = rng.getRandomInt(0,_exposedPositions.seq.size()-1);
    string mutAA = possibleMutations.substr(mutIndex,1);
    string previousAA = current.seq.substr(seqPosIndex,1);
    current.seq.replace(seqPosIndex,1,mutAA);

    // Eval
    saltBridgeResult sbResult = scoreSaltBridgePositions(current,_sbScoreTable,_sys);
    MSLOUT.stream() << "SCORE: "<<sbResult.score<<" wt: "<<wt_sbresult.score<<endl;

    if (sbResult.score == 0.0){
      current.seq = wtSeq;
      continue;
    }

    // MC test for Acceptance:
    if (MCMngr.accept(sbResult.score)){

      //cout << "New SEQ: "<<current.seq<<" score: "<<sbResult.score<<endl;
      //cout << "WT  SEQ: "<<wtSeq<<" score: "<<wt_sbresult.score<<endl;
      bool keepIt = false;

      map<string, bool>::iterator it;

      it = mcDataMap.find(sbResult.sbId);

      if (it == mcDataMap.end()){

	mcDataMap[sbResult.sbId] = true;

	if (_mcData.size() >= _opt.numSequenceModels){

	  if (sbResult.score < _mcData.top().first.score){

	    // Remove highest score, then add
	    _mcData.pop();
	    _mcData.push(pair<saltBridgeResult,subSeq>(sbResult,current));
	    keepIt = true;
	  }

	} else {

	  _mcData.push(pair<saltBridgeResult,subSeq>(sbResult,current));
	  keepIt = true;
	}

      }


      if (!keepIt){
	current.seq.replace(seqPosIndex,1,previousAA);
      }
    } else {

      current.seq.replace(seqPosIndex,1,previousAA);
   }
  }
  MSLOUT.stream() << "MC Completed -> "<<MCMngr.getReasonCompleted()<<endl;

}

void runModelSaltBridges(
			 System &_sys,
			 Options &_opt,
			 std::priority_queue< std::pair<saltBridgeResult,subSeq>, std::vector< std::pair<saltBridgeResult,subSeq> >, compareScores> &_mcData){


  // TODO: Use sorted rotamer library..

  System baseModel = _sys;

  ofstream fastaFile;
  if (_opt.fasta_only){
    //cout << "Opening: "<<MslTools::stringf("RSB_%s.fasta", MslTools::getFileName(_opt.pdb).c_str())<<endl;
    fastaFile.open(MslTools::stringf("RSB_%s.fasta", MslTools::getFileName(_opt.pdb).c_str()));
  }
    
  // Loop over solutions
  uint model = 0;
  vector<pair<saltBridgeResult,subSeq> > sbResultInOrder;
  while (!_mcData.empty() && model <= _opt.numStructuralModels){
    sbResultInOrder.push_back(_mcData.top());
    _mcData.pop();
    model++;
  }

  
  model = 0;
  for (int sb = sbResultInOrder.size()-1; sb >= 0; sb--){

    //fprintf(stdout, "WT : %s %8.3f\n", wtSeq.c_str(),wt_sbresult.score);
    //fprintf(stdout, "SEQ: %s %8.3f\n",mcData.top().second.seq.c_str(),mcData.top().first.score);
    saltBridgeResult sbr = sbResultInOrder[sb].first;

    MSLOUT.fprintf(stdout, "SOLUTION: %8.3f\n",sbr.score);

    if (_opt.fasta_only){
      vector<pair<string,string> > posIds;
      for (uint i = 0; i < sbr.mutations.size();i++){
	//cout << "PosID1: "<<getPosIdAndMutation(baseModel, sbr.mutations[i].first,_opt).first<<endl;
	//cout << "PosID2: "<<getPosIdAndMutation(baseModel, sbr.mutations[i].second,_opt).first<<endl;
	posIds.push_back(getPosIdAndMutation(baseModel, sbr.mutations[i].first,_opt));
	posIds.push_back(getPosIdAndMutation(baseModel, sbr.mutations[i].second,_opt));
      }
      //cout << "Writing in fasta file\n";
      fastaFile << MslTools::stringf(">%s%sRSB_%04d\n", MslTools::getFileName(_opt.pdb).c_str(),_opt.appendName.c_str(),model);
      for (uint c = 0; c < baseModel.chainSize();c++){
	for (uint r = 0; r < baseModel.getChain(c).positionSize();r++){
	  Residue &res = baseModel.getChain(c).getResidue(r);
	  string aa = MslTools::getOneLetterCode(res.getResidueName());
	  for (uint i = 0; i < posIds.size();i++){
	    if (res.getPositionId() == posIds[i].first){
	      aa = MslTools::getOneLetterCode(posIds[i].second);
	      break;
	    }
	  }
	  if (aa != "X")
	    fastaFile << MslTools::stringf("%1s",aa.c_str());
	  }
	fastaFile << endl;
      }
      model++;
      continue;
    }
    
    vector<int> variablePositions;
    for (uint i = 0; i < sbr.mutations.size();i++){
      MSLOUT.fprintf(stdout,"\t%s %s\n",sbr.mutations[i].first.c_str(),sbr.mutations[i].second.c_str());

      // Mutate System, return position id of mutated position , blank string if no mutation...
      string posId = mutate(baseModel, sbr.mutations[i].first,_opt);
      if (posId != ""){
	MSLOUT.fprintf(stdout,"\t\t Mutating1 %s to %s\n", posId.c_str(),sbr.mutations[i].first.c_str());
	variablePositions.push_back(baseModel.getPositionIndex(posId));
      }
      
      posId = mutate(baseModel, sbr.mutations[i].second,_opt);
      if (posId != ""){
	MSLOUT.fprintf(stdout,"\t\t Mutating2 %s to %s\n", posId.c_str(),sbr.mutations[i].second.c_str());
	variablePositions.push_back(baseModel.getPositionIndex(posId));
      }
    }

    // Change HIS to HSD
    for (uint i = 0 ; i < baseModel.positionSize();i++){
      if (baseModel.getPosition(i).getResidueName() == "HIS"){
	baseModel.getPosition(i).setResidueName("HSD");
      }
    }

//    // Fasta output
//    if (_opt.fasta_only){
//
//      fastaFile << MslTools::stringf(">RSB_%04d_%s\n",model, MslTools::getFileName(_opt.pdb).c_str());
//      for (uint c = 0; c < baseModel.chainSize();c++){
//	for (uint r = 0; r < baseModel.getChain(c).positionSize();r++){
//	  Residue &res = baseModel.getChain(c).getResidue(r);
//	  string aa = MslTools::getOneLetterCode(res.getResidueName());
//	  if (aa != "X")
//	    fastaFile << MslTools::stringf("%1s",aa.c_str());
//	  }
//	fastaFile << endl;
//      }
//      model++;
//      baseModel = _sys;
//      continue;
//    }
    
    MSLOUT.stream() << " Quenching..."<<endl;
    // Quench model..
    Quench quencher(_opt.topfile, _opt.parfile, _opt.rotlib);    
    quencher.setVariableNumberRotamers(50,10); // Big residues, Small residues
    System quenchedSys = quencher.runQuench(baseModel,variablePositions);

    // Write out model..
    quenchedSys.writePdb(MslTools::stringf("%s_resurfaced_%06d.pdb", MslTools::getFileName(_opt.pdb).c_str(),model+1));

    // Increment the models
    model++;

    baseModel = _sys;
  }

  if (_opt.fasta_only){
    //cout << "Closing fasta file\n";
    fastaFile.close();
  }
  
  MSLOUT.stream() << "Done salt-bridge modeling"<<endl;

}

pair<string,string> getPosIdAndMutation(System &_sys, string _mutationId, Options &_opt){
  string chainId;
  int resnum;
  string icode;
  string old_identity;
  string new_identity;
  if (!MslTools::parseMutationId(_mutationId,chainId,resnum,icode,old_identity,new_identity)){
    cerr << "ERROR 1234 Couldn't parse mutation: "<<_mutationId<<endl;
    exit(1234);
  }

  return pair<string,string>(MslTools::getPositionId(chainId, resnum,icode),new_identity);
}
string mutate(System &_sys, string _mutationId, Options &_opt){      

  string chainId;
  int resnum;
  string icode;
  string old_identity;
  string new_identity;
  if (!MslTools::parseMutationId(_mutationId,chainId,resnum,icode,old_identity,new_identity)){
    cerr << "ERROR 1234 Couldn't parse mutation: "<<_mutationId<<endl;
    exit(1234);
  }
  if (old_identity == new_identity){
    return "";
  }
  if (old_identity == "HIS") old_identity = "HSD";
  if (new_identity == "HIS") new_identity = "HSD";

  string posId = MslTools::getPositionId(chainId, resnum,icode);
 
  PDBTopology pdbTop;
  pdbTop.readRotamerLibrary(_opt.rotlib);
  pdbTop.setAddAtomsFromRotLib(true);

  // Get backbone atoms  
  AtomPointerVector backboneAtoms = pdbTop.getBackboneAtoms(_sys.getPosition(posId).getCurrentIdentity());

  // Mutate
  AtomContainer newAtoms = pdbTop.getResidue(MslTools::getIdentityId(chainId,resnum,icode,new_identity),backboneAtoms,1);
  _sys.getPosition(posId).addIdentity(newAtoms.getAtomPointers(),new_identity);
  _sys.getPosition(posId).setActiveIdentity(new_identity);
  _sys.getPosition(posId).removeIdentity(old_identity);


  return posId;
}
