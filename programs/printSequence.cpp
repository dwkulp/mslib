/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2012 The MSL Developer Group (see README.TXT)
 MSL Libraries: http://msl-libraries.org

If used in a scientific publication, please cite: 
 Kulp DW, Subramaniam S, Donald JE, Hannigan BT, Mueller BK, Grigoryan G and 
 Senes A "Structural informatics, modeling and design with a open source 
 Molecular Software Library (MSL)" (2012) J. Comput. Chem, 33, 1645-61 
 DOI: 10.1002/jcc.22968

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
#include "MslTools.h"
#include "OptionParser.h"
#include "release.h"

#include "System.h"
#include "Residue.h"

using namespace std;
using namespace MSL;

#include "printSequence.h"

int main(int argc, char *argv[]) {

	Options opt = setupOptions(argc, argv);

	System sys;
	sys.readStructureFile(opt.pdb);// could be PDB or CIF
	if (opt.fasta)
	  fprintf(stdout, ">%s\n",MslTools::getFileName(opt.pdb).c_str());
	for (uint c = 0; c < sys.chainSize();c++){
	  for (uint r = 0; r < sys.getChain(c).positionSize();r++){
	      Residue &res = sys.getChain(c).getResidue(r);
	      string aa = MslTools::getOneLetterCode(res.getResidueName());
	      if (aa != "X")
		fprintf(stdout,"%1s",aa.c_str());
	  }
      fprintf(stdout,"\n");
	}

}


Options setupOptions(int theArgc, char * theArgv[]){
	// Create the options
	Options opt;

	// Parse the options
	OptionParser OP;
	OP.readArgv(theArgc, theArgv);
	OP.setRequired(opt.required);	
	OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option


	if (OP.countOptions() == 0){
		cout << "Usage: printSequence " << endl;
		cout << endl;
		cout << "\n";
		cout << "pdb PDB\n";
		cout << endl;
		exit(0);
	}

	opt.pdb = OP.getString("pdb");
	if (OP.fail()){
		cerr << "ERROR 1111 no pdb specified."<<endl;
		exit(1111);
	}
	opt.fasta = OP.getBool("fasta");

	return opt;
}
