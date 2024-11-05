/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2013 The MSL Developer Group (see README.TXT)
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

#include "AtomSelection.h"
#include "System.h"
#include "OptionParser.h"
#include "release.h"
#include "AtomSelection.h"
#include "Transforms.h"
#include "RegEx.h"
#include "MslExceptions.h"
#include "MslTools.h"

using namespace std;

using namespace MSL;


string programName = "alignMolecules";
string programDescription = "This programs aligns two PDB based on a subset of atoms";
string programAuthor = "Alessandro Senes";
string programVersion = "1.0.4";
string programDate = "24 April 2013";
string mslVersion =  MSLVERSION;
string mslDate = MSLDATE;

 // variables to store time at start and end of application
time_t start_time, finish_time;

/******************************************
 *      =======  STRUCTURE Options  =======
 * 
 *  This structure contains all the user
 *  defined options
 * 
 *  The function
 *     Options parseOptions(int _argc, char * _argv[], Options defaults)
 *  parses the argument (using the OptionParser object) and return
 *  a filled Options structure (opt) which is used during the program
 *
 *  Another Options structure (defaults) stores some default
 *  values that will be given to opt if the user did not input them as an
 *  argument
 * 
 ******************************************/
struct Options {
	string commandName; // name of this program

	/***** OPTION VARIABLES START HERE ******/
	// MOLECULE
	string pdb1; // first pdb
	string pdb2; // second pdb
	vector<string> sele1; // atom selection one
	vector<string> sele2; // atom selection two
        vector<string> regex; // Regular expresssion to match primary sequence
	unsigned int model1; // model of pdb 1 (if NMR structure)
	bool setModel1; // if model1 was given activates a setter
	unsigned int model2; // model of pdb 2 (if NMR structure)
	bool setModel2; // if model1 was given activates a setter
	string outputPdb2; // the new pdb
	bool writeAllModels; // determes if all models are writted for a multi model input PDB
	bool noAlign; // only calculate the rmsd between the selections, do not move the molecules
	bool noOutputPdb; // do not write an output pdb
	string outputdir;  // the directory with the output for the run
	string configfile;  // name of the configuration file

	/***** MANAGEMENT VARIABLES ******/
	string pwd; // the present working directory obtained with a getenv
	string host; // the host name obtained with a getenv
	bool version; // ask for the program version
	bool help; // ask for the program help

	ofstream * cout_fs; // file stream for the redirected standard output
	ofstream * cerr_fs; // file stream for the redirected standard error
	streambuf * coutbuf; // pointer for storing the cout streambuffer
	streambuf * cerrbuf; // pointer for storing the cerr streambuffer
	streambuf * coutsbuf; // pointer to the filestr streambuffer
	streambuf * cerrsbuf; // pointer to the filestr streambuffer
	bool coutRedirected; // true if cout was redirected

	bool errorFlag; // true if there are errors
	bool warningFlag; // true if there are warnings
	string errorMessages; // error messages
	string warningMessages; //warning messages

	vector<string> allowed; //list of allowed options
	vector<string> required; //list of required options

	vector<string> disallowed;  // disallowed options that were given
	vector<string> missing; // required options that were not given
	vector<string> ambiguous; // required options that were not given
	vector<string> defaultArgs; // the default arguments can be specified in command line without "--option"
	vector<vector<string> > equivalent; // this links short options to long ones (for example -x can be given for --extended)

	string OPerrors; //the errors from the option parser



};




/******************************************
 *  
 *  =======  PREDECLARATION OF THE FUNCTIONS  =======
 *
 ******************************************/
Options parseOptions(int _argc, char * _argv[], Options defaults);
void usage();
void version();
void help(Options defaults);
void addAtomsByRegEx(AtomPointerVector &_input, AtomPointerVector &_output,string &_regex);



/******************************************
 *  
 *         =======  MAIN  =======
 *
 ******************************************/
int main(int argc, char *argv[]) {
	
	// store start time
	time(&start_time);

	cout << "Start: " << ctime(&start_time);
	cout << "Program: " << programName << endl;
	cout << "Program description: " << programDescription << endl;
	cout << "Program author: " << programAuthor << endl;
	cout << "Program version: " << programVersion << " " << programDate << endl;
	cout << "MSL version: " << mslVersion << " " << mslDate << endl;
	/******************************************************************************
	 *                          === SETTINGS THE DEFAULTS ===
	 *
	 *  Put here the defaults for some options that are
	 *  not always required
	 ******************************************************************************/
	Options defaults;
	defaults.noAlign = false;
	defaults.noOutputPdb = false;
	defaults.writeAllModels = false;
	

	/******************************************************************************
	 *                             === OPTION PARSING ===
	 *
	 *  Parse the command line (and possibly input file)
	 *  options. It will also create an input file based
	 *  on the current options in the output directory
	 *  that can be used to re-run the program with the same
	 *  configuration
	 ******************************************************************************/
	Options opt = parseOptions(argc, argv, defaults);
	if (opt.errorFlag) {
		cerr << endl;
		cerr << "The program terminated with errors:" << endl;
		cerr << endl;
		cerr << opt.errorMessages << endl;
		cerr << endl;
		cerr << opt.OPerrors << endl;

		//printErrors(opt);
		usage();
		exit(1);
	}


	cout << "Read pdb 1: " << opt.pdb1 << endl;
	System sys1;
	if (!sys1.readPdb(opt.pdb1)) {
		cerr << "Unable to open pdb " << opt.pdb1 << endl;
		exit(1);
	}
	if (opt.setModel1) {
		// set the current NMR model
		sys1.setActiveModel(opt.model1);
	}

	cout << "Read pdb 2: " << opt.pdb2 << endl;
	System sys2;
	if (!sys2.readPdb(opt.pdb2)) {
		cerr << "Unable to open pdb " << opt.pdb2 << endl;
		exit(1);
	}
	if (opt.setModel2) {
		// set the current NMR model
		sys2.setActiveModel(opt.model2);
	}


	AtomPointerVector av1 = sys1.getAtomPointers();
	AtomPointerVector av2 = sys2.getAtomPointers();

	AtomPointerVector alignAtoms1;
	if (opt.sele1.size() == 0) {

	  // select all atoms
	  alignAtoms1.insert(alignAtoms1.end(), av1.begin(), av1.end());

	} else {
		AtomSelection sel1(av1);
		for (unsigned int i=0; i<opt.sele1.size(); i++) {
			char c [1000];
			sprintf(c, "keyatoms, %s", opt.sele1[i].c_str());
			AtomPointerVector selAtom = sel1.select(c);
			alignAtoms1.insert(alignAtoms1.end(), selAtom.begin(), selAtom.end());
		}
	}

	int numAts_priorTo_RegEx = alignAtoms1.size();
	if (opt.regex.size() != 0){
	  AtomPointerVector newRegExAts1;

	  for (uint r = 0; r < opt.regex.size();r++){
	    cout << "REGEX["<<r<<"]: ->"<< opt.regex[r]<<"<-"<<endl;
	    addAtomsByRegEx(alignAtoms1,newRegExAts1,opt.regex[r]);
	  }
	  
	  alignAtoms1.clear();
	  alignAtoms1 = newRegExAts1;	  
	}

	cout << "->Selected " << alignAtoms1.size() << " (prior to RegEx: "<<numAts_priorTo_RegEx<<") reference atoms for pdb " << opt.pdb1 << endl;


	AtomPointerVector alignAtoms2;
	if (opt.sele2.size() == 0) {

	  // select all atoms
	  alignAtoms2.insert(alignAtoms2.end(), av2.begin(), av2.end());

	} else {
		AtomSelection sel2(av2);
		for (unsigned int i=0; i<opt.sele2.size(); i++) {
			char c [1000];
			sprintf(c, "keyatoms, %s", opt.sele2[i].c_str());
			AtomPointerVector selAtom = sel2.select(c);
			alignAtoms2.insert(alignAtoms2.end(), selAtom.begin(), selAtom.end());
		}
	}

	/*

	numAts_priorTo_RegEx = alignAtoms2.size();
	if (opt.regex != ""){
	  AtomPointerVector newRegExAts2;
	  addAtomsByRegEx(alignAtoms2,newRegExAts2,opt);
	  alignAtoms2.clear();
	  alignAtoms2 = newRegExAts2;
	}
	*/
	numAts_priorTo_RegEx = alignAtoms2.size();
	if (opt.regex.size() != 0){
	  AtomPointerVector newRegExAts2;

	  for (uint r = 0; r < opt.regex.size();r++){
	    addAtomsByRegEx(alignAtoms2,newRegExAts2,opt.regex[r]);
	  }
	  
	  alignAtoms2.clear();
	  alignAtoms2 = newRegExAts2;	  
	}

	
	cout << "->Selected " << alignAtoms2.size() << " (prior to RegEx: "<<numAts_priorTo_RegEx<<") reference atoms for pdb " << opt.pdb2 << endl;

	if (alignAtoms1.size() != alignAtoms2.size()) {
		cerr << "The number of atoms selected for pdb 1 (" << alignAtoms1.size() << ") does not match the number of atoms selected for pdb 2 (" << alignAtoms2.size() << ")" << endl;
		exit(1);
	}

	cout << endl;
	cout << "Set 1 ================================" << endl;
	cout << endl;
	cout << alignAtoms1;
	cout << endl;
	cout << "Set 2 ================================" << endl;
	cout << endl;
	cout << alignAtoms2;
	cout << endl;
	cout << "      ================================" << endl;

	if (opt.noAlign) {
		// only calc the rmsd
		double rmsd = alignAtoms2.rmsd(alignAtoms1);
		cout << endl;
		cout << opt.pdb2 << " was NOT aligned to " << opt.pdb1 << "." << endl;
		cout << "RMSD " << rmsd << endl;
	} else {
		Transforms tm;
		if (!tm.rmsdAlignment(alignAtoms2,alignAtoms1,av2)) {
			cerr << "Alignment failed!" << endl;
			exit(1);
		}
		//double rmsd = tm.getRMSD();
		double rmsd = alignAtoms2.rmsd(alignAtoms1);
		Matrix rotMatrix = tm.getLastRotationMatrix();
		CartesianPoint translation = tm.getLastTranslation();
		cout << rotMatrix << endl;
		cout << translation << endl;
		cout << endl;
		cout << "Aligned " << opt.pdb2;
		if (opt.setModel2) {
			cout << " (model " << opt.model2 << ")";
		}
		cout << " to " << opt.pdb1;
		if (opt.setModel1) {
			cout << " (model " << opt.model1 << ")";
		}
		cout << "." << endl;
		cout << "RMSD " << rmsd << endl;
		if (!opt.noOutputPdb) {
			if (!sys2.writePdb(opt.outputPdb2, opt.writeAllModels)) {
				cerr << "Cannot open " << opt.outputPdb2 << " for writing" << endl;
				exit(1);
			}
			cout << "Written to " << opt.outputPdb2 << endl;
		} else {
			cout << "Aligned PDB not written" << endl;
		}
	}



	return 0;

}


void addAtomsByRegEx(AtomPointerVector &_input, AtomPointerVector &_output,string &_regex){

  // Only use 1 atom per residue, "CA"
  AtomSelection sel(_input);
  AtomPointerVector CAats = sel.select("name CA");

  RegEx re;
  re.setStringType(RegEx::PrimarySequence); 
  vector<pair<int,int> > matchingResidueIndices = re.getResidueRanges(CAats,_regex);

  // Do we get first, last or other match?
	    //for (uint r = matchingResidueIndices[matchingResidueIndices.size()-1].first; r <= matchingResidueIndices[matchingResidueIndices.size()-1].second;r++){
  for (uint i = 0; i < matchingResidueIndices.size();i++){
      for (uint r = matchingResidueIndices[i].first; r <= matchingResidueIndices[i].second;r++){
	Atom *at = CAats[r];
	_output.push_back(at);
      }
  }

}

Options parseOptions(int _argc, char * _argv[], Options defaults) {

	/******************************************
	 *  Pass the array of argument and the name of
	 *  a configuration file to the ArgumentParser
	 *  object.  Then ask for the value of the argument
	 *  and collect error and warnings.
	 *
	 *  This function returns a Options structure
	 *  defined at the head of this file 
	 ******************************************/
	
	Options opt;

	/******************************************
	 *  Set the allowed and required options:
	 *
	 *  Example of configuartion file:
	 *  
	 ******************************************/
	vector<string> required;
	vector<string> allowed;

	opt.required.push_back("pdb1");
	opt.required.push_back("pdb2");
	opt.allowed.push_back("sele1");
	opt.allowed.push_back("sele2");
	opt.allowed.push_back("regex");
	opt.allowed.push_back("model1");
	opt.allowed.push_back("model2");
	opt.allowed.push_back("outputPdb2");
	opt.allowed.push_back("writeAllModels");
	opt.allowed.push_back("noAlign");
	opt.allowed.push_back("noOutputPdb");
	opt.allowed.push_back("outputdir");
	opt.allowed.push_back("configfile");
	opt.allowed.push_back("version"); // --version
	opt.allowed.push_back("v"); // -v is equivalent to --version
	opt.allowed.push_back("help"); // --help
	opt.allowed.push_back("h"); // -h is equivalent to --help

	opt.equivalent.push_back(vector<string>());
	opt.equivalent.back().push_back("v");
	opt.equivalent.back().push_back("version");
	opt.equivalent.push_back(vector<string>());
	opt.equivalent.back().push_back("h");
	opt.equivalent.back().push_back("help");

	opt.defaultArgs.push_back("configfile");


	//OptionParser OP(_argc, _argv);
	OptionParser OP;
	OP.setShortOptionEquivalent(opt.equivalent);
	OP.readArgv(_argc, _argv);
	OP.setDefaultArguments(opt.defaultArgs); // a pdb file value can be given as a default argument without the --pdbfile option
	OP.setRequired(opt.required);
	OP.setAllowed(opt.allowed);
	OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"
	if (OP.countOptions() == 0) {
		usage();
		exit(0);
	}
	opt.configfile = OP.getString("configfile");
	if (opt.configfile != "") {
		OP.readFile(opt.configfile);
		if (OP.fail()) {
			opt.errorFlag = true;
			opt.errorMessages += "Cannot read configuration file " + opt.configfile + "\n";
			exit(1);
		}
	}

	/*****************************************
	 *  VERSION AND HELP
	 *
	 *  --version or -v arguments print the version number
	 *  --help prints the usage and help
	 *****************************************/
	opt.version = OP.getBool("version");
	//if (OP.fail()) {
	//	opt.version = OP.getBool("v");
	//}

	if (opt.version) {
		version();
		exit(0);
	}

	opt.help = OP.getBool("help");
//	if (OP.fail()) {
//		opt.help = OP.getBool("h");
//	}

	if (opt.help) {
		help(defaults);
		exit(0);
	}

	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	if (!OP.checkOptions()) {
		opt.errorFlag = true;
		opt.OPerrors = OP.getErrors();
		return opt;
	}
	opt.errorFlag = false;
	opt.warningFlag = false;


	/*****************************************
	 *  CHECK THE GIVEN OPTIONS
	 *****************************************/
	//opt.commandName = OP.getCommandName();
	//if (!OP.checkOptions()) {
	//	opt.errorFlag = true;
	//	opt.disallowed = OP.getDisallowedOptions();
	//	opt.missing = OP.getMissingOptions();
	//	opt.ambiguous = OP.getAmbiguousOptions();
	//	return opt;
	//}
	opt.errorFlag = false;
	opt.warningFlag = false;


	opt.pdb1 = OP.getString("pdb1");
	if (OP.fail()) {
		opt.errorMessages = "Option name of first pdb file \"pdb1\" not specified";
		opt.errorFlag = true;
	}

	opt.pdb2 = OP.getString("pdb2");
	if (OP.fail()) {
		opt.errorMessages = "Option name of second pdb file \"pdb2\" not specified";
		opt.errorFlag = true;
	}

	opt.outputPdb2 = OP.getString("outputPdb2");
	if (OP.fail()) {
		opt.warningMessages = "Option name of output aligned second pdb file \"outputPdb2\" not specified";
		opt.warningFlag = true;
		string base = MslTools::pathTail(opt.pdb2);
		base = MslTools::pathRoot(base);
		opt.outputPdb2 = base + (string)"-aligned.pdb";
	}

	opt.writeAllModels = OP.getBool("writeAllModels");
	if (OP.fail()) {
		opt.writeAllModels = defaults.writeAllModels;
	}

	int index = 0;
	while (true) {
		string sele = OP.getString("sele1", index);
		if (OP.fail()) {
			break;
		}
		opt.sele1.push_back(sele);
		index++;
	}
	//if (opt.sele1.size() == 0) {
	//	opt.sele1.push_back("all");
	//}

	index = 0;
	while (true) {
		string sele = OP.getString("sele2", index);
		if (OP.fail()) {
			break;
		}
		opt.sele2.push_back(sele);
		index++;
	}
	//if (opt.sele2.size() == 0) {
	//	opt.sele2.push_back("all");
	//}
	index = 0;
	while (true) {
		string sele = OP.getString("regex", index);
		if (OP.fail()) {
			break;
		}
		opt.regex.push_back(sele);
		index++;
	}

	opt.model1 = OP.getInt("model1");
	if (!OP.fail()) {
		opt.setModel1 = true;
	}

	opt.model2 = OP.getInt("model2");
	if (!OP.fail()) {
		opt.setModel2 = true;
	}

	opt.noAlign = OP.getBool("noAlign");
	if (OP.fail()) {
		opt.noAlign = defaults.noAlign;
	}
	opt.noOutputPdb = OP.getBool("noOutputPdb");
	if (OP.fail()) {
		opt.noOutputPdb = defaults.noOutputPdb;
	}
	// return the Options structure
	return opt;

}

void usage() {
	cout << endl;
	cout << "Run as" << endl;
	cout << "   % " << programName << " --configfile <config.txt>" << endl;
	cout << "For help" << endl;
	cout << "   % " << programName << " -h" << endl;
	cout << endl;
}

void version() {
	cout << endl;
	cout << "Program " << programName << " v." << programVersion << ", " << programDate << ", (MSL v." << mslVersion << " " << mslDate << ")" << endl;
}

void help(Options defaults) {
       	cout << "Run  as:" << endl;
	cout << " % alignMolecules --pdb1 <target.pdb> [--model1 <N>] --pdb2 <movable.pdb> [--model2 <N>] [--sele1 <atom selection> --sele2 <atom selection> [--sele1 <atom selection> --sele2 <atom selection>]] [--outputPdb2 <output.pdb>] [--noAlign]" << endl;
	cout << endl;
	cout << " *******************************************************************************************" << endl;
	cout << " * Options:                                                                                *" << endl;
	cout << " *     pdb1             : name of the pdb to align to (target)                             *" << endl;
	cout << " *     pdb2             : name of the pdb to be aligned (movable)                          *" << endl;
	cout << " *     model1           : optional model number for NMR multi-model PDBs (default #1)      *" << endl;
	cout << " *     model2           : optional model number for NMR multi-model PDBs (default #1)      *" << endl;
	cout << " *     writeAllModels   : if PDB 2 is multi-model, write all models (default is false)     *" << endl;
	cout << " *     sele1            : selection (one or more) of atoms from PDB 1 to be used for the   *" << endl;
	cout << " *                        alignment                                                        *" << endl;
	cout << " *     sele2            : selection (one or more) of atoms from PDB 2 to be used for the   *" << endl;
	cout << " *                        alignment (the total number of atoms must match)                 *" << endl;
	cout << " *     noAlign          : calculate the current RMSD without aligning                      *" << endl;
	cout << " *     noOutputPdb      : do not write a PDB file out                                      *" << endl;
	cout << " *                                                                                         *" << endl;
	cout << " *     regex            : regular expression(s) for selections (order agnostic)            *" << endl;
	cout << " * NOTE: No order is assumed WITHIN a selection but if multiple --sele1/sele2 are          *" << endl;
	cout << " *       given order is preserved.                                                         *" << endl;
	cout << " *       If no selection is given, all atoms are used.                                     *" << endl;
	cout << " ******************************************************************************************" << endl;
	cout << endl;
}

