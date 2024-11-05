/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)n
 Copyright (C) 2009 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan

This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, 
 USA, or go to http://www.gnu.org/copyleft/lesser.txt.
----------------------------------------------------------------------------
*/
#include <string>
#include <map>
#include <fstream>
#include <signal.h>
#include "OptionParser.h"
#include "Timer.h"
#include "System.h"
#include "AtomSelection.h"
#include "getMeasurement.h"

using namespace std;

using namespace MSL;

// Global objects
Timer t;
double startTime = 0.0;

int main(int argc, char *argv[]){

	// Option Parser
	Options opt = setupOptions(argc,argv);

	// Create a system from the structural input options
	startTime = t.getWallTime();

	//read in list of PDBs to compare to first PDB
	vector<string> list;
	if (MslTools::pathExtension(opt.pdb) == "pdb"){
	  list.push_back(opt.pdb);
	} else {
	  MslTools::readTextFile(list,opt.pdb);
	}
    
	for (uint i = 0; i < list.size();i++){

	  System sys;
	  sys.readPdb(list[i]);

	  if (opt.atom1 != "" && opt.atom2 != ""){
	    double dist = sys.getAtom(opt.atom1).distance(sys.getAtom(opt.atom2));
	    fprintf(stdout, "%40s %10s %10s %8.3f\n",list[i].c_str(),opt.atom1.c_str(),opt.atom2.c_str(),dist);
	  }

	  if (opt.atom1sel != "" && opt.atom2sel != ""){
	      AtomSelection sel(sys.getAtomPointers());
	      AtomPointerVector ats1 = sel.select(opt.atom1sel);
	      AtomPointerVector ats2 = sel.select(opt.atom2sel);
	      
	      for (uint a1 = 0; a1 < ats1.size();a1++){
		for (uint a2 = 0; a2 < ats2.size();a2++){
		  double dist = ats1(a1).distance(ats2(a2));
		  fprintf(stdout, "%40s %10s %10s %8.3f\n",list[i].c_str(),ats1(a1).getAtomId().c_str(),ats2(a2).getAtomId().c_str(),dist);
		}
	      }

	  }

	  if (opt.terminiDist){

	    PDBReader *pin = sys.getPDBReader();
	    map<string,double> boundingCoords = pin->getBoundingCoordinates();
	    double xcoor = boundingCoords["maxX"] - boundingCoords["minX"];
	    double ycoor = boundingCoords["maxY"] - boundingCoords["minY"];
	    double zcoor = boundingCoords["maxZ"] - boundingCoords["minZ"];

	    int numAAchains = 0;
	    for (uint c = 0; c < sys.chainSize();c++){
	      if (sys.getChain(c).positionSize() > 20){
		if (sys.getChain(c).getPosition(0).atomExists("CA") || sys.getChain(c).getPosition((int(sys.getChain(c).positionSize()/2))).atomExists("CA") || sys.getChain(c).getPosition(sys.getChain(c).positionSize()-1).atomExists("CA")) {
		  numAAchains++;
		}
	      }
	    }
	    for (uint c = 0; c < sys.chainSize();c++){
	      AtomSelection sel(sys.getChain(c).getAtomPointers());
	      AtomPointerVector ats = sel.select("name CA");
	      if (ats.size() > 1){
		double dist = ats(0).distance(ats(ats.size()-1));
		cout << MslTools::stringf("%50s,%1s,%6.2f,%8d, %6.3f,%6.0f, %6.0f, %6.0f, %8d\n",list[i].c_str(),sys.getChain(c).getChainId().c_str(),dist,sys.getChain(c).positionSize(),dist/(double)sys.getChain(c).positionSize(),xcoor,ycoor,zcoor,numAAchains);
	      }
	    }
	    
	  }
	  

	  


	}

}




Options setupOptions(int theArgc, char * theArgv[]){
	Options opt;

	OptionParser OP;


	OP.setRequired(opt.required);
	OP.setAllowed(opt.optional);
	OP.readArgv(theArgc, theArgv);

	if (OP.countOptions() == 0){
		cout << "Usage:" << endl;
		cout << endl;
		cout << "getMeasurement --pdb PDB"<<endl;
		exit(0);
	}
	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 pdb not specified.\n";
		exit(1111);
	}
	opt.atom1 = OP.getString("atom1");
	opt.atom2 = OP.getString("atom2");
	opt.atom1sel = OP.getString("atom1sel");
	opt.atom2sel = OP.getString("atom2sel");
	opt.terminiDist = OP.getBool("terminiDistance");
	
	return opt;
}







