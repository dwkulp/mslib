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
#ifndef BACKRUB_H
#define BACKRUB_H

#include <vector>
#include <string>
#include "Chain.h"
#include "AtomPointerVector.h"
#include "System.h"
#include "RandomNumberGenerator.h"

namespace MSL { 
class BackRub {

	public:
		BackRub();
		~BackRub();
		
		void localSample(Chain &_ch, int _startResIndex, int _endResIndex, int _numFragments);
		std::vector<AtomContainer *> multiSample(Chain &_ch, int _startResIndex, int _endResIndex,  int _numSamplingAttempts, int _numModels, double _rmsdCutoff=2.0);
		bool doRub(Chain &_ch, int startRes, int endRes);
		
		AtomPointerVector&  getAtomPointers() { return sys.getAllAtomPointers(); }
		std::string getNMRString() { return sysNMRFormat; }
	private:
		
		void doMinorRotation(Residue &_r1, Residue &_r2, CartesianPoint &_targetAtom1, CartesianPoint &_targetAtom2);

		System sys;
		RandomNumberGenerator rng;
		std::string sysNMRFormat;
};

}

#endif
