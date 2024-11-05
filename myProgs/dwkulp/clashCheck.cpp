
#include "MslTools.h"
#include "OptionParser.h"
#include "release.h"
#include "PhiPsiStatistics.h"
#include "clashCheck.h"
#include "System.h"
#include "Chain.h"
#include "Residue.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;
using namespace MSL;
using namespace MslTools;

int main(int argc, char *argv[]) {

    Options opt = setupOptions(argc, argv);

    vector<string> lines;
    if (MslTools::pathExtension(opt.pdblist).compare("pdb") == 0){
      cout << "JUST A PDB!\n";
      lines.push_back(opt.pdblist);
    } else {
      cout << "PDBLIST IS: ."<<MslTools::pathExtension(opt.pdblist)<<"."<<endl;
      MslTools::readTextFile(lines,opt.pdblist);
    }


    for (uint i = 0; i < lines.size();i++){


	System sys;
	sys.readPdb(lines[i]);

        map<string,int> clashMap;
        map<string,vector<string> > clashingPositions;
        int total_clashes = 0;

	// Each chain
	for (uint c1 = 0; c1 < sys.chainSize();c1++){

	  for (uint p1 = 0; p1 < sys.getChain(c1).positionSize();p1++){

	    // Only 
	    bool skipIt = true;
	    if (opt.residueTypes.size() > 0) {
	      for (uint rt = 0; rt< opt.residueTypes.size();rt++){
		if (sys.getChain(c1).getPosition(p1).getResidueName() == opt.residueTypes[rt]){
		  skipIt = false;
		}
	      }
	    } else {
	      skipIt = false;
	    }

	    if (skipIt) {
	      continue;
	    }

	    int pos_clashes = 0;
	    for (uint c2 = 0; c2 < sys.chainSize();c2++){
	      for (uint p2 = 0; p2 < sys.getChain(c2).positionSize();p2++){
		if (abs ((int)p1 - (int)p2) < 3) continue;

		bool skipIt2 = true;
		if (opt.residueTypes.size() > 0) {
		  for (uint rt = 0; rt< opt.residueTypes.size();rt++){
		    if (sys.getChain(c2).getPosition(p2).getResidueName() == opt.residueTypes[rt]){
		      skipIt2 = false;
		    }
		  }
		} else {
		  skipIt2 = false;
		}

		if (skipIt2) {
		  continue;
		}

	      
		for (uint a1 = 0; a1 < sys.getChain(c1).getPosition(p1).atomSize();a1++){
		  Atom &at1 = sys.getChain(c1).getPosition(p1).getAtom(a1);
		  if (at1.getElement() == "H") continue;
		  for (uint a2 = a1+1; a2 < sys.getChain(c2).getPosition(p2).atomSize();a2++){
		    Atom &at2 = sys.getChain(c2).getPosition(p2).getAtom(a2);
		    if (at2.getElement() == "H") continue;

		    double dist = at1.distance(at2);

		    if (dist < opt.dist){
		    //cout << "POSITIONS ARE CLASHING: "<<sys.getChain(c1).getPosition(p1).getPositionId()<<" "<<sys.getChain(c1).getPosition(p1).getResidueName()<<" - "<<sys.getChain(c2).getPosition(p2).getPositionId()<<" "<<sys.getChain(c2).getPosition(p2).getResidueName()<<endl<<std::flush;
		      pos_clashes++;
		      clashingPositions[sys.getChain(c1).getPosition(p1).getPositionId()].push_back(sys.getChain(c2).getPosition(p2).getPositionId());
		    }
		  
		  
		  }
		}
	      }
	    }
	    //	    cout << "POS_CLASH: "<<pos_clashes<<endl<<std::flush;
	    if (pos_clashes > 0){
	      clashMap[sys.getChain(c1).getPosition(p1).getPositionId()] = pos_clashes;
	    }
	    total_clashes += pos_clashes;

	  }


	  
	}

	if (opt.printClashes){
	  map<string,int>::iterator it;
	  for (it = clashMap.begin(); it != clashMap.end();it++){
	    fprintf(stdout, "%10s %8d",it->first.c_str(),it->second);
	    for (uint cp=0; cp < clashingPositions[it->first].size();cp++){
	      fprintf(stdout, ",%10s", clashingPositions[it->first][cp].c_str());
	    }
	    fprintf(stdout,"\n");
	  }
	}

	fprintf(stdout, "TOTAL: %s %8d\n",lines[i].c_str(),total_clashes);

	if (total_clashes < opt.tooManyClashes) {
	  fprintf(stdout,"MIN CLASHES: %s",lines[i].c_str());
	}

    }


}

Options setupOptions(int theArgc, char * theArgv[]){
    // Create the options
    Options opt;

    // Parse the options
    OptionParser OP;
    OP.setRequired(opt.required);	
    OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option
    OP.readArgv(theArgc, theArgv);

    if (OP.countOptions() == 0){
	cout << "Usage: clashCheck " << endl;
	cout << endl;
	cout << "\n";
	cout << "pdblist PDB\n";
	cout << endl;
	exit(0);
    }

    opt.pdblist = OP.getString("pdblist");
    if (OP.fail()){
	cerr << "ERROR 1111 no pdblist specified."<<endl;
	exit(1111);
    }
    opt.tooManyClashes = OP.getInt("tooManyClashes");
    if (OP.fail()){
	opt.tooManyClashes = 0;
    }

    opt.residueTypes = OP.getStringVector("resTypes");
    if (OP.fail()){
    }
    opt.dist  = OP.getDouble("dist");
    if (OP.fail()){
      opt.dist = 1.5;
    }
    opt.printClashes = OP.getBool("reportClashes");
    
    return opt;
}
