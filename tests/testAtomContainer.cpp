/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
 Sabareesh Subramaniam, Ben Mueller

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
#include "AtomContainer.h"
#include "Transforms.h"
#include "MslTools.h"

using namespace MSL;
using namespace std;

int main(int argc, char *argv[]) {

	// the program requires the location of the "exampleFiles" as an argument
	if (argc < 2) {
		cout << "USAGE:\nexample_AtomContainer_usage <path_of_exampleFiles_directory> <atom identifier>\n\nAtom identifier: either A,37,CB or A,37,ILE,CB (with insertion code A,37A,CB)" << endl;
		exit(0);
	}

	AtomContainer molecule;
	molecule.readPdb(argv[1]);
	cout << molecule;

	if (argc < 3) {
		cout << endl;
		cout << "MISSING ATOM IDENTIFIER!\n\nUSAGE:\nexample_AtomContainer_usage <path_of_exampleFiles_directory> <atom identifier>\n\nAtom identifier: either A,37,CB or A,37,ILE,CB (with insertion code A,37A,CB)" << endl;
		exit(0);
	} else {
		if (molecule.atomExists(argv[2])) {
			cout << "Atom " << argv[2] << " exists" << endl;
			cout << "Print it using getLastFoundAtom()" << endl;
			cout << "   " << molecule.getLastFoundAtom() << endl;
			cout << "Print it using molecule(" << argv[2] << ")" << endl;
			cout << "   " << molecule(argv[2]) << endl;
		} else {
			cout << "Atom " << argv[2] << " does NOT exists" << endl;
		}
	}

	return 0;
}
