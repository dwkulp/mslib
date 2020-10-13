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


#ifndef CIFWRITER_H
#define CIFWRITER_H
/*
  This class handles writing of PDB files.
 */

// MSL Includes
#include "Writer.h"
#include "MslTools.h"
#include "CIFFormat.h"

// Storage formats
#include "CartesianPoint.h"
#include "AtomPointerVector.h"
#include "FormatConverter.h"

// STL Includes
#include <vector>



namespace MSL { 

class CIFWriter : public Writer {

	public:
		// Constructors/Destructors
		CIFWriter();
		CIFWriter(const std::string &_filename);
		virtual ~CIFWriter();

		// Get/Set

		// Member Functions
		bool write(AtomPointerVector &_av, bool _noHydrogens=false, bool _addHeader=true, bool _addTail=true);

		bool open();               // There is a default implementation
		bool open(const std::string &_filename); // There is a default implementation
		bool open(const std::string &_filename, int mode); // There is a default implementation
		bool open(std::stringstream &_ss);
		void close();


	protected:		
	private:

};

//Inlines go HERE
inline bool CIFWriter::open() {bool success = Writer::open(); return success;}
inline bool CIFWriter::open(const std::string &_filename) {bool success = Writer::open(_filename); return success;}
inline bool CIFWriter::open(const std::string &_filename, int mode) {bool success = Writer::open(_filename, mode); return success;}
inline bool CIFWriter::open(std::stringstream &_ss) {fileHandler = stringstyle; bool success = Writer::open(_ss); return success;}


}

#endif