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

#ifndef CIFFORMAT_H
#define CIFFORMAT_H

//MSL Includes
#include "MslTools.h"
#include "Atom.h"
#include "AtomPointerVector.h"

//STL Includes
#include <string>
#include <iostream>
#include <cstring>
#include <sstream>

namespace MSL { 
class CIFFormat {
  
/*
  CIF Formating information goes here.  As well as a simple structure for containing a single line
  of a CIF.
*/
       public:
                static Atom * createAtomFromAtomSiteLine(const std::string &_pdbAtomLine, std::vector<std::string> &_fields, bool _allowPartialRead=false,bool _useAuthFields=true);

                static std::string createLoopAtomSite(AtomPointerVector &_ats, bool _addHeader=true, bool _addTail=true);
                static std::string createLoopAtomSite(AtomPointerVector &_ats, std::vector<std::string> &_fields, bool _addHeader=true, bool _addTail=true);

                static std::string createAtomSiteLineFromAtom(Atom &_at,unsigned int _atomNumber=0);
                static std::string createAtomSiteLineFromAtom(Atom &_at, std::vector<std::string> &_fields,unsigned int _atomNumber=0);


	protected:
	private:


};


}

#endif
