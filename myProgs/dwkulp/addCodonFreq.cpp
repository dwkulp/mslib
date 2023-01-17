#include <iostream>
#include <cstdlib>
//#include <queue>
//#include <algorithm>
#include <map>
#include <string>

#include "System.h"
#include "OptionParser.h"
#include "MslOut.h"
#include "SysEnv.h"
#include "PDBWriter.h"
#include "FastaReader.h"

#include "addCodonFreq.h"

using namespace std;
using namespace MSL;


// MslOut 
static MslOut MSLOUT("addCodonFreq");
static SysEnv SYSENV;

  
int main(int argc, char *argv[]) {

    // Parse commandline options
  Options opt = setupOptions(argc,argv);

  // Read in original pdb
  System sys;
  sys.readStructureFile(opt.pdb);

  // opt.sequence is the DNA seq corresponding to a structure file

  if (opt.fasta != ""){
    FastaReader fin;
    fin.open(opt.fasta);
    fin.read();
    fin.close();
    map<string,string> seqs = fin.getSequences();
    for (map<string,string>::iterator it = seqs.begin(); it != seqs.end();it++){
      opt.seq = it->second;
    }
  }
  // Freq table should be Codon, AA, Freq, Freq.thousand (AGC, Ser, 0.55, 0.42)
  vector<string> lines;
  MslTools::readTextFile(lines, opt.table);

  // Map codon to freq
  map<string, double> codonFreq;
  map<string, string> codon2AA;
  for (uint i = 0; i < lines.size();i++){
    if (lines[i].substr(0,1) == "#") continue;
    if (MslTools::trim(lines[i]) == "") continue;
    vector<string> toks = MslTools::tokenize(lines[i]," ");
    //cout << "Tokens: "<<toks[0]<<","<<toks[1]<<","<<toks[2]<<endl;
    
    codonFreq.insert(pair<string,double>(toks[0],std::stod(toks[2].c_str())));
    codon2AA.insert(pair<string,string>(toks[0],MslTools::toUpper(toks[1])));

  }
  cout << "Now add codon freq..\n";
  int seqIndex = 0;
  for (uint c = 0; c < sys.chainSize();c++){
    Chain &ch = sys.getChain(c);
    
    int startIndex = 0;
    int endIndex = ch.positionSize();

    for (uint r = startIndex; r < endIndex;r++){
      Residue &res   = ch.getResidue(r);
      string aa      = res.getResidueName();
      string codon   = opt.seq.substr(seqIndex, 3);
      string codonAA = codon2AA[codon];
      //cout << "OUT: sequence index: "<<seqIndex<<" codon: "<<codon<<" codonAA: "<<codonAA<<" != "<<res.toString()<<endl;
      
      // Double check that we have the same AA.
      if (aa == codonAA){
	res.getAtomPointers().setTempFactor(codonFreq[codon]);
      } else {
	cerr << "ERROR: sequence index: "<<seqIndex<<" codon: "<<codon<<" codonAA: "<<codonAA<<" != "<<res.toString()<<endl;
	exit(1234);
      }

      // Increment to next codon
      seqIndex += 3;
    }
    
  }

  sys.writePdb(opt.out);
  
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
    cout << "addCodonFreq --pdb pdb --seq DNASEQ --table foo.txt  --out foo.pdb\n";

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

  opt.seq = OP.getString("seq");
  if (OP.fail()){
    opt.seq="";
  }
  opt.fasta = OP.getString("fasta");
  if (OP.fail()){
    opt.fasta="";
  }

  if (opt.seq == "" && opt.fasta == ""){
    cerr << "ERROR 1111 seq and fasta not specified.\n";
    exit(1111);
  }
  opt.table = OP.getString("table");
  if (OP.fail()){
    string dir = "MSL_DIR";
    stringstream ss;
    if (SYSENV.isDefined(dir)){
	ss << SYSENV.getEnv("MSL_DIR") << "/tables/mammalian.txt";
	if (!MslTools::fileExists(ss.str())){
	  cerr << "ERROR 1111 default MSL Env table does not exist ("<<ss.str()<<").\n";
	  exit(1111);
	}
    } else {
      cerr << "ERROR 1111 table and default MSL Env table not specified.\n";
      exit(1111);
    }
    opt.table =	ss.str();
    cerr << "Table set to default: "<<ss.str()<<endl;
  }

  opt.out = OP.getString("out");
  if (OP.fail()){
    opt.out = "tmp.pdb";
  }
  
  return opt;
}
