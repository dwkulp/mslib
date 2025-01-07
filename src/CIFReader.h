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

#ifndef CIFREADER_H
#define CIFREADER_H
/*

 */

//MSL Includes
#include "Reader.h"
#include "MslTools.h"
#include "CIFFormat.h"

// Storage Formats
#include "CartesianPoint.h"
#include "AtomPointerVector.h"


// STL Includes
#include <vector>

/**
 * This class will provide an object which is able
 * to read in and interpret CIF files.
 */
namespace MSL { 
class CIFReader : public Reader {

	public:
		// Constructors/Destructors
		CIFReader();
		CIFReader(const std::string &_filename);
        CIFReader(const CIFReader & _reader);
		CIFReader(std::stringstream &_stream);
		~CIFReader();

		bool read(bool _noHydrogens=false, bool _allowPartialRead=false);
		bool read(std::string &_inputString);

		AtomPointerVector & getAtomPointers(); 
		unsigned int size() const;
		Atom * operator[](unsigned int _n);


		void reset();

		unsigned int getNumberOfModels() const;

		// Two fields of data one is 'label' (set by PDB) and the other is 'auth' (set by author), use author by default
		void setUseAuthFields(bool _useAuthFields) {useAuthFields = _useAuthFields;}

	protected:		
	private:
		void deletePointers();

		void processData(std::string &_category, std::vector<std::string> &_fields, std::string &_line);
		AtomPointerVector atoms;

		bool singleAltLocFlag;

		std::vector<Matrix *> symmetryRotations;
		std::vector<CartesianPoint *> symmetryTranslations;

		std::vector<Matrix *> biounitRotations;
		std::vector<CartesianPoint *> biounitTranslations;

		Matrix *scaleRotation;
		CartesianPoint *scaleTranslation;
		
		std::vector<double> unitCellParams;

		std::map<std::string,double> boundingCoords;

		unsigned int numberOfModels;
		bool useAuthFields;

};

//Inlines go HERE
/**
* This method will return a std::vector of atoms found in this CIF file.
*
* @return A std::vector or atoms from the PDB file.
*/
inline AtomPointerVector& CIFReader::getAtomPointers() { return atoms; }

/**
 * This method will delete all data held in the CIFReader.  All
 * Atom pointers that were previously saved off will no longer be valid
 * after this reset.
 */
inline void CIFReader::reset() {deletePointers();}



/*
  Why do I need this even though I have one in Reader.h
 */
inline bool CIFReader::read(std::string &_inputString){
	fileName = "string";
	stringStreamPtr = new std::stringstream();
	stringStreamPtr->str(_inputString);
	fileHandler = stringstyle;

	// Now call read..
	return read();
}

inline unsigned int CIFReader::getNumberOfModels() const {return numberOfModels;}


}

#endif
