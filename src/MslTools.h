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


#ifndef MSLTOOLS_H_
#define MSLTOOLS_H_

// STL includes
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <map>
#include <limits>
#include <cstdlib>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <time.h>
#include <stdio.h>
#include <stdarg.h>

// MSL includes
#include "MslExceptions.h"
#include "Real.h"
#include "release.h"
#include "RandomNumberGenerator.h"

// BOOST Includes
#ifdef __BOOST__
#include <boost/regex.hpp>
#endif



namespace MSL{
    namespace MslTools {




	/*
             ******************************************
	     *          CONSTANTS
	     ******************************************
	*/

	const float  floatMax  = std::numeric_limits<float>::max();
	const double doubleMax = std::numeric_limits<double>::max();
	const int    intMax    = std::numeric_limits<int>::max();

	// Molar Gas Constant in kcal/mol/K -  Based on "CODATA recommended values of the fundamental physical constants:2010"
	// Reviews of modern physics, vol 84, oct-dec 2012
	const double R    = 1.972e-3;

	


	/*
             ******************************************
	     *          STRING FUNCTIONS
	     ******************************************
	*/

	// trim trailing and leading whitespaces
	std::string trim(const std::string & _str, const std::string & _trimString=" \t\n\r");
	std::vector<std::string> trim(const std::vector<std::string> & _str, const std::string & _trimString=" \t\n\r");

	// Replace
	bool replace(std::string &_string, const std::string &_replace, const std::string &_with, bool _replaceAll=false);

	// CONVERSIONS
	double toDouble(const std::string & _string, const std::string & _msg=(std::string)"");
        Real toReal(const std::string & _string, const std::string & _msg=(std::string)"");
	int toInt(const std::string & _string, const std::string & _msg=(std::string)"");
	unsigned int toUnsignedInt(const std::string & _string, const std::string & _msg=(std::string)"");
	bool toBool(const std::string & _string, const std::string & _msg=(std::string)"");
	std::string intToString(const int & _i, const std::string & _msg=(std::string)"");
	std::string unsignedIntToString(const unsigned int & _i, const std::string & _msg=(std::string)"");
	std::string doubleToString(const double & _d, const std::string & _msg=(std::string)"");
	std::string stringf(const char * _format, ...);

	// CHECK CONTENT OF STRING
	bool isDigitChars(std::string _input); // 0-9
	bool isAlphaNumericChars(std::string _input); // 0-9 A-Z a-z
	bool isAlphaChars(std::string _input); // A-Z a-z
	bool isWhiteSpaces(std::string _input); // " " \t \n \r

	// String comparisions
	unsigned int hamming_distance(const std::string &_str1, const std::string &_str2);


	 // split a string made by digits followed by letters into an int and a string (i.e. "734B" -> 734 and "B")
	bool splitIntAndString(const std::string & _input, int & _intResult, std::string & _stringResult);

	// ATOM AND RESIDUE IDs
	// The Atom Id is in the form of "A 37 CA" or "A 37A CA" with an insertion code
	std::string getAtomId(std::string _chainid, int _resnum, std::string _icode, std::string _atomName, unsigned int _skiplevel=0);
	bool parseAtomId(std::string _atomId, std::string & _chainid, int & _resnum, std::string & _icode, std::string & _atomName, unsigned int _skiplevels=0);
	bool compareAtomIds(std::string _id1, std::string _id2, unsigned int _skiplevels=0);

	// The Residue Id is in the form of "A 37" or "A 37A" with an insertion code
	std::string getPositionId(std::string _chainid, int _resnum, std::string _icode, unsigned int _skiplevel=0);
	//std::string getPositionId(std::string _chainid, int _resnum, unsigned int skiplevel=0); // blank icode
	bool parsePositionId(std::string _posId, std::string & _chainid, int & _resnum, std::string & _icode, unsigned int _skiplevels=0);
	bool comparePositionIds(std::string _id1, std::string _id2, unsigned int _skiplevels=0);

	// The Identity Id is in the form of "A 37 ILE" or "A 37A ILE" with an insertion code
	std::string getIdentityId(std::string _chainid, int _resnum, std::string _icode, std::string _identity, unsigned int _skiplevel=0);
	//std::string getIdentityId(std::string _chainid, int _resnum, std::string _identity, unsigned int skiplevel=0); // blank icode
	bool parseIdentityId(std::string _residueId, std::string & _chainid, int & _resnum, std::string & _icode, std::string & _identity, unsigned int _skiplevels=0);
	bool compareIdentityIds(std::string _id1, std::string _id2, unsigned int _skiplevels=0);

	// The Rotamer Id is in the form of "A 37 ILE 3" or "A 37A ILE 3" with an insertion code
	std::string getRotamerId(std::string _chainid, int _resnum, std::string _icode, std::string _identity, unsigned int _conformation);
	bool parseRotamerId(std::string _rotamerId, std::string & _chainid, int & _resnum, std::string & _icode, std::string & _identity, unsigned int &_conformation);

	// The Atom of Identity Id is in the form of "A 37 ILE CA" or "A 37A ILE CA" with an insertion code
	std::string getAtomOfIdentityId(std::string _chainid, int _resnum, std::string _icode, std::string _identity, std::string _atomName, unsigned int _skiplevel=0);
	bool parseAtomOfIdentityId(std::string _atomId, std::string & _chainid, int & _resnum, std::string & _icode, std::string & _identity, std::string & _atomName, unsigned int _skiplevels=0);
	bool compareAtomOfIdentityIds(std::string _id1, std::string _id2, unsigned int _skiplevels=0);

	// Identity Id in the form "A 37 ILE ARG" or "A,37,ILE,ARG"
	bool parseMutationId(std::string _mutationId, std::string & _chainid, int & _resnum, std::string & _icode, std::string & _identity,std::string &_newIdentity, unsigned int _skiplevels=0);

	// case
	std::string toUpper(const std::string & _input);
	std::string toLower(const std::string & _input);

	// JOIN AND SPLIT
	std::vector<std::string> tokenize(const std::string & _input, const std::string & _delimiter=" ", bool _allowEmtpy=false);
	std::vector<std::vector<std::string> > dualLevelTokenizer(const std::string & _input, const std::string & _delimiter=" ", const std::string & _leftBraket="[", const std::string _rightBraket="]");
	std::vector<std::string> extractBraketed(const std::string & _input, const std::string & _leftBraket="[", const std::string _rightBraket="]");
	//std::vector<std::string> tokenizeWord(const std::string & _input, const std::string & _delimiter=" ", const bool & _trim=false);
	std::vector<std::string> tokenizeAndTrim(const std::string & _input, const std::string & _delimiter=" ", bool _allowEmtpy=false, const std::string & _trimString=" \t\n\r");

	std::string joinLines(const std::vector<std::string> & _input, const std::string & _spacer=" ");

	std::vector<std::string> joinConnectedLines(const std::vector<std::string> & _input, std::string _marker="\\", std::string _spacer=" ");
	//std::vector<std::string> joinBackslashedLines(const std::vector<std::string> & _input, const std::string & _spacer=" ");

	std::vector<std::string> removeEmptyLines(const std::vector<std::string> & _input);
	std::string uncomment(const std::string & _input, const std::string & _commentString="#");
	std::vector<std::string> uncomment(const std::vector<std::string> & _input, const std::string & _commentString="#");
	

	// Extract the name of a file from a full path , subtract extension as well.
	bool readTextFile(std::vector<std::string> & _container, const std::string & _filename);

	bool fileExists(std::string _filename);
	bool directoryExists(std::string _filename);

	// Parsing file without chain ids (e.g. 1ABC )
	bool getPdbsFromFile(std::string _filename, std::vector<std::string> &_pdbs, std::string _pdbPath, std::string _cifPath);

	// Parsing file with PDBIDCHIAINID format (e.g. 1ABCA ) note chain id CAN be more than one character, this code will handle that.
	bool getPdbsFromFile(std::string _filename, std::vector<std::string> &_pdbs, std::vector<std::string> &_chainIds, std::string _pdbPath, std::string _cifPath, bool _parseChain=true);

	std::string getMSLversion();

	// RegEx Functions  (only works with compile __BOOST__ flag on, otherwise returns false immediately)
	//   Remember escaped characters need to be double-escaped in expression variable so '\s' is '\\s' 
	bool regex(std::string lineToMatch, std::string expression, std::vector<std::string> &matches);
	bool regex(std::string lineToMatch, std::string expression, std::vector<std::pair<int,int> > &matches);
		
	/*
             ******************************************
	     *          MATH FUNCTIONS
	     ******************************************
	*/

	double correlate(std::vector<double> _data1, std::vector<double> _data2);
	double smartRound(double coord, double gridSize);
	double mod(double x, double y);
	double setPrecision(double _d, unsigned int _significantDigits);

	Real round(Real value);

	/*
             ******************************************
	     *          BOLTZAMM PROBABILITIES
	     ******************************************
	*/
	// the energy of an ensemble of states E = sum p(i) * E(i) for all states i
	double getBoltzmannEnsembleEnergy(double _temp, std::vector<double>& _energies);
	std::vector<double> getBoltzmannProbabilities(double _temp, std::vector<double>& _energies);

	
	/*
             ******************************************
	     *          PATH FUNCTIONS
	     ******************************************
	*/
	std::string pathRoot(std::string _path); // /the/path/theFile.ext -> /the/path/theFile
	std::string pathHead(std::string _path); // /the/path/theFile.ext -> /the/path
	std::string pathExtension(std::string _path); // /the/path/theFile.ext -> ext
	std::string pathTail(std::string _path); // /the/path/theFile.ext -> theFile.ext
	std::string getFileName(std::string fullpath); // /the/path/theFile.ext -> theFile

	std::string createDir(std::string _name);
	std::string outputFileNameParser(std::string _name);
	bool   mkNestedDir(std::string _dir, mode_t _mode); // WHAT IS MODE?
	std::string getRandomAlphaNumString(unsigned int _size, bool _alphaOnly=false);
	//unsigned int getRandomInt(unsigned int _max);

	/*
              ******************************************
              * Tools to convert different AA codes.
	      ******************************************
	*/
	static std::map<std::string, std::string> threeToOneLetter;
	static std::map<std::string, std::string> oneToThreeLetter;
	void loadAAConversionTables();
	std::string getOneLetterCode(std::string threeLetterCode);
	std::string getThreeLetterCode(std::string oneLetterCode);


	/*
              ******************************************
              * STL Helper functions (sort predicates,etc..)
	      ******************************************
        */

	bool sortPairIntDoubleAscending(const std::pair<int,double> &left, const std::pair<int,double> &right);
	bool sortPairIntDoubleDecending(const std::pair<int,double> &left, const std::pair<int,double> &right);

	/*
              ******************************************
              * Help sorting residues by resnum and icode
	      ******************************************
        */
	bool sortByResnumIcodeAscending(int _resnum1, std::string _icode1, int _resnum2, std::string _icode2);

        template <class A, class B> inline bool mapHasKey(std::map<A, B> &_mapObj, A _key) {
            return( _mapObj.find(_key) != _mapObj.end());
        };


	/*
              ******************************************
              * Sorting a std::vector of doubles
	      ******************************************
	*/
	void quickSort(std::vector<double> & _vec);
	void quickSort(std::vector<double> & _vec, unsigned int _left, unsigned int _right);
	double partition(std::vector<double> & _vec, unsigned int _left, unsigned int _right, unsigned int _pivotIndex);
	/*
              ******************************************
              * Sorting a std::vector of doubles and keeping track of
	      * the original index in another std::vector (_index
	      * should initially be _index[i] := i)
	      ******************************************
	*/
	void quickSortWithIndex(std::vector<double> & _vec, std::vector<unsigned int> &_index);
	double partitionWithIndex(std::vector<double> & _vec, unsigned int _left, unsigned int _right, unsigned int _pivotIndex, std::vector<unsigned int> &_index);
	void quickSortWithIndex(std::vector<double> & _vec, unsigned int _left, unsigned int _right, std::vector<unsigned int> &_index);

	/*
             ******************************************
	     *         VECTOR<DOUBLE> FUNCTIONS 
	     ******************************************
	*/
	void normalizeVector(std::vector<double> & _vec); /*!< \brief Normalize the vector so that the sum is 1 */
	void normalizeCumulativeVector(std::vector<double> & _vec); /*!< \brief Normalize so that the last element is 1 */


	/*
             ******************************************
	     *          COLOR FUNCTIONS
	     ******************************************
	*/
	std::vector<double> getRGB(std::vector<double> &_startRGB, std::vector <double> &_endRGB, double _minValue, double _maxValue, double _value);
	void rgb2hsv(std::vector<double> &_rgb, std::vector<double> &_hsv);
	void hsv2rgb(std::vector<double> &_hsv, std::vector<double> &_rgb);


	std::string getBaseMSLDirectory();

    }
};


// inlined functions

#endif
