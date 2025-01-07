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


#include "CIFReader.h"

using namespace MSL;
using namespace std;

#include "MslOut.h"
static MslOut MSLOUT("CIFReader");

/**
 * Simple constructor.
 */
CIFReader::CIFReader() : Reader() {
  	numberOfModels = 0;
	singleAltLocFlag = false;
	scaleTranslation = new CartesianPoint(0.0,0.0,0.0);
	scaleRotation = new Matrix(3,3,0.0);
	(*scaleRotation)[0][0] = 1.0;
	(*scaleRotation)[1][1] = 1.0;
	(*scaleRotation)[2][2] = 1.0;
	useAuthFields = true;
}
/**
 * With this constructor the user specifies the filename
 * of the PDB to be read.
 *
 * @param _filename  The name of the PDB file to be read.
 */
CIFReader::CIFReader(const std::string &_filename) : Reader(_filename) {
	numberOfModels = 0;
	singleAltLocFlag = false;
	scaleTranslation = new CartesianPoint(0.0,0.0,0.0);
	scaleRotation = new Matrix(3,3,0.0);
	(*scaleRotation)[0][0] = 1.0;
	(*scaleRotation)[1][1] = 1.0;
	(*scaleRotation)[2][2] = 1.0;
	useAuthFields = true;
}
/**
 * A copy constructor.  All of the atoms from the given CIFReader are
 * copied into the new CIFReader.
 *
 * @param _reader The CIFReader to be copied.
 */
CIFReader::CIFReader(const CIFReader & _reader) {
	for (AtomPointerVector::const_iterator k=_reader.atoms.begin(); k!= _reader.atoms.end(); k++) {
		atoms.push_back(new Atom(**k));
	}
	numberOfModels = _reader.numberOfModels;
	singleAltLocFlag = _reader.singleAltLocFlag;
	scaleTranslation = new CartesianPoint(0.0,0.0,0.0); 
	scaleRotation = new Matrix(3,3,0.0);
	useAuthFields = true;
}
/**
 * A constructor which will read input data from a std::stringstream.
 *
 * @param _ss The std::stringstream to get data from.
 */
CIFReader::CIFReader(std::stringstream &_ss) : Reader(_ss)     {
	numberOfModels = 0;
	singleAltLocFlag = false;
	scaleTranslation = new CartesianPoint(0.0,0.0,0.0); 
	scaleRotation = new Matrix(3,3,0.0);
	useAuthFields = true;
}

/**
 * The deconstructor.  All data will be deleted, so any Atom pointers
 * that were previously saved off will no longer be valid after the CIFReader
 * object has been destroyed.
 */
CIFReader::~CIFReader() {
	deletePointers();
	close();
}


/**
 * This method will delete all of the Atom pointers
 * held in this CIFReader.
 */
void CIFReader::deletePointers() {
	for (AtomPointerVector::iterator k=atoms.begin(); k!=atoms.end(); k++) {
		delete *k;
	}

	atoms.clear();

}

/**
 * This method will actually read the data in from the
 * given CIF file.  Currently the reader only looks at 
 * _atom_site information.  All other information 
 * simply ignored.
 */
bool CIFReader::read(bool _noHydrogens, bool _allowPartialRead) {

	if (!is_open()) {
		return false;
	}

	MSLOUT.stream() << "FILE IS OPEN\n";
	try { 


		numberOfModels = 0;
		int inLoop = 0;
		string inLoop_category = "";
		vector<string> inLoop_fields;
		
		while (!endOfFileTest()){
			string line = Reader::getLine();

			
			// check the length
			string header = "";

			// Blank line, also indicates end of loop_
			if (line.substr(0,1) == "#"){

			  //MSLOUT.stream() << "COMMENT LINE:"<<line<<endl;

			  inLoop_category = "";
			  inLoop_fields.clear();
			  inLoop = 0;
			  
			  continue;
			}
			
			// start of loop_ section
			if (line.substr(0,5) == "loop_"){
			  inLoop = 1;
			  //MSLOUT.stream() << "In a _loop" <<endl;
			  continue;
			}

			// Find category first line in loop
			if (inLoop == 1 && line.substr(0,1) == "_"){
			  
			  vector<string> toks = MslTools::tokenize(line,".");
			  inLoop_category = MslTools::trim(toks[0]);

			  inLoop_fields.push_back(MslTools::trim(toks[1]));

			  inLoop = 2;
			  //MSLOUT.stream() << "First Line of _loop: " <<toks[1]<<endl;
			  continue;
			}

			if (inLoop == 2 && line.substr(0,1) == "_"){

			  vector<string> toks = MslTools::tokenize(line,".");
			  inLoop_fields.push_back(MslTools::trim(toks[1]));
			  //MSLOUT.stream() << "Next Line of _loop: " <<toks[1]<<endl;
			  continue;
			}

			// Data line in _loop
			if (inLoop == 2 && line.substr(0,1) != "_"){
			  //MSLOUT.stream() << "Data in _loop" <<endl;
			  processData(inLoop_category, inLoop_fields, line);
			  continue;
			}
			

		}
		
	} catch(...){
		cerr << "ERROR 5623 in CIFReader::read()\n";
		return false;
	}

	return true;
}



void CIFReader::processData(string &_category, vector<string> &_fields, string &_line){

  if (_category == "_atom_site"){
    MSLOUT.stream() << "CREATE ATOM BUDDY!:"<<_line<<endl;
	Atom *a = CIFFormat::createAtomFromAtomSiteLine(_line, _fields, useAuthFields);
    atoms.push_back(a);
  }
  

}
