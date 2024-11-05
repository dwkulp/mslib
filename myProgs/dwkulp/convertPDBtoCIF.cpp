#include <iostream>
#include <cstdlib>
#include <queue>
#include <algorithm>
#include <map>

#include "System.h"
#include "OptionParser.h"
#include "MslOut.h"
#include "SysEnv.h"
#include "PDBWriter.h"

#include "convertPDBtoCIF.h"

using namespace std;
using namespace MSL;


// MslOut 
static MslOut MSLOUT("convertPDBtoCIF");
static SysEnv SYSENV;

  
int main(int argc, char *argv[]) {

    // Parse commandline options
  Options opt = setupOptions(argc,argv);

  // Read in original pdb
  System sys;
  sys.readPdb(opt.pdb);

  CIFWriter cifout;
  cifout.open(MslTools::stringf("%s.cif",MslTools::getFileName(opt.pdb).c_str()));

  vector<string> chainIds = {"A","B","C","D","E","F","G","H","J","I","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"};
  int numChains = sys.chainSize();
  int numModels = sys.getNumberOfModels();
  for (uint i = 0; i < numModels;i++){
    cout << "Working on model: "<<i<<endl;
    
    string newChainAd = "";
    double fract = i / chainIds.size();
    double result = floor(fract);
    int repeat = ((int)result)+1;    
    int index = i % chainIds.size();
    
    char char_array[2];
    strcpy(char_array, chainIds[index].c_str());
    newChainAd = std::string(repeat, char_array[0]);
    sys.setActiveModel(i);
    AtomContainer ac(sys.getAtomPointers());
    System newSys(ac.getAtomPointers());

    for (uint c = 0; c < numChains;c++){
      string newChain = MslTools::stringf("%s%s",newChainAd.c_str(),newSys.getChain(c).getChainId().c_str());
      cout << "\tModel "<<i<<" "<<newSys.getChain(c).getChainId()<<" "<<newChain<<endl;
      newSys.getChain(c).setChainId(newChain);
    }

    if (i == 0){
      cifout.write(newSys.getAtomPointers(),false,true,false);
    } else {
      if (i != numModels-1){
	cifout.write(newSys.getAtomPointers(),false,false,false);
      } else {
	cifout.write(newSys.getAtomPointers(),false,false,true);
      }
    }
  }

  cifout.close();
  
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
    cout << "simpleNanoparticleModeling --pdb pdb\n";

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
