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
#include <numeric>    //inner_product
#include <functional> //plus, equal_to, not2
#include <cstdlib>
#include <iostream>

using namespace MSL;
using namespace std;



string MslTools::trim(const string & _str, const string &_trimString){
	string out = _str;

	// trim leading whitespace
	string::size_type  notwhite = out.find_first_not_of(_trimString);
	out.erase(0,notwhite);

	// trim trailing whitespace
	notwhite = out.find_last_not_of(_trimString); 
	out.erase(notwhite+1); 

	return out;
}

vector<string> MslTools::trim(const vector<string> & _str,const string &_trimString){
	vector<string> out;
	for (vector<string>::const_iterator k=_str.begin(); k!=_str.end(); k++) {
		out.push_back(trim(*k,_trimString));
	}
	return out;
}

double MslTools::toDouble(const string & _string, const string & _msg) {
	istringstream ss(_string);
	double d = 0.0;
	ss >> d;

	if (!ss) throw ConvertDoubleException(" string "+ _string +" is not a double. " + _msg);

	return d;
}

Real MslTools::toReal(const string & _string, const string & _msg) {
	istringstream ss(_string);
	Real d = 0.0;
	ss >> d;

	if (!ss) throw ConvertDoubleException(" string "+ _string +" is not a Real. " + _msg);

	return d;
}

unsigned int MslTools::toUnsignedInt(const string & _string, const string & _msg) {

	istringstream ss(_string);
	unsigned int i = 0;
	ss >> i;

	if (!ss) throw ConvertIntException(" string "+ _string +" is not an unsigned int. "+ _msg);

	return i;
}

int MslTools::toInt(const string & _string, const string & _msg) {

	istringstream ss(_string);
	int i = 0;
	ss >> i;

	if (!ss) throw ConvertIntException(" string "+ _string +" is not an int. "+ _msg);

	return i;
}

bool MslTools::toBool(const string & _string, const string & _msg) {
	string upperStr= _string;
	upperStr= MslTools::toUpper(upperStr);
	if (upperStr== "0" || upperStr== "FALSE" || upperStr== "F") {
		// default is true unless 0, false or f are given
		return false;
	}
	return true;
}

string MslTools::intToString(const int & _i, const string & _msg){
	stringstream ss;
	ss << _i;

	return ss.str();
	
}
string MslTools::unsignedIntToString(const unsigned int & _i, const string & _msg){
	stringstream ss;
	ss << _i;

	return ss.str();
	
}
string MslTools::doubleToString(const double & _d, const string & _msg){
	stringstream ss;
	ss << _d;

	return ss.str();
	
}

string MslTools::toUpper(const string & _input) {
	string out = _input;
	for(int i=0; i<out.length(); i++) {
		out[i] = toupper(out[i]);
	}
	return out;
}

string MslTools::toLower(const string & _input) {
	string out = _input;
	for(int i=0; i<out.length(); i++) {
		out[i] = tolower(out[i]);
	}
	return out;
}

vector<string> MslTools::tokenize(const string & _input, const string & _delimiter, bool _allowEmtpy){
	vector<string> results;

	if (_input == "") {
		if (_allowEmtpy) {
			results.push_back(_input);
		}
		return results;
	}
	
	if (_allowEmtpy) {
		size_t prePos = 0;
		size_t pos  = _input.find(_delimiter);
		unsigned int delimiterSize = _delimiter.size();
		string left = _input, right;

		while (pos != std::string::npos) {
			results.push_back(left.substr(prePos, pos));
			if( pos + delimiterSize <= left.size() ) {
				left = left.substr(pos + delimiterSize, left.size() );
			} else {
				left = "";
			}
			pos  = left.find(_delimiter);
		}

		results.push_back(left);
	} else {
		int start  = _input.find_first_not_of(_delimiter);
		int end    = 0;
		string cur = _input;


		while (start != std::string::npos){
			end    = _input.find_first_of(_delimiter, start);
			results.push_back(_input.substr(start, end-start));
			start  = _input.find_first_not_of(_delimiter, end);
		}
	}
	return results;

}

vector<string> MslTools::tokenizeAndTrim(const string & _input, const string & _delimiter, bool _allowEmtpy, const string & _trimString) {
	vector<string> results = tokenize(_input, _delimiter, _allowEmtpy);
	return trim(results, _trimString);
}

/*
vector<string> MslTools::tokenizeWord(const string & _input, const string & _delimiter, const bool & _trim){
	vector<string> results;

	int start  = _input.find_first_not_of(_delimiter);
	int end    = 0;
	int loop   = 0;
	while ( start != std::string::npos){
		if (loop == 0) start = 0;
		end    = _input.find(_delimiter, start);

		if (_trim){
			results.push_back(trim(_input.substr(start, end-start)));
		} else {
			results.push_back(_input.substr(start, end-start));			
		}
		start  = _input.find_first_not_of(_delimiter,end);
		loop++;
	}

	return results;
}
*/


vector<vector<string> > MslTools::dualLevelTokenizer(const string & _input, const string & _delimiter, const string & _leftBraket, const string _rightBraket) {
	/*********************************************************************
	 *  This function turns a string like this
	 *        A1 B1 C1 [D1 D2 D3] [E1 E2] F1 [G1 G2]
	 *  into a vector of vector
	 *    [0][0] A1
	 *    [1][0] B1
	 *    [2][0] C1
	 *    [3][0] D1
	 *    [3][1] D2
	 *    [3][2] D3
	 *    [4][0] E1
	 *    [4][1] E2
	 *    [5][0] F1
	 *    [6][0] G1
	 *    [6][1] G2
	 *  
	 *  The bracket type can be specified (defaul is currently []
	 *********************************************************************/


	vector<vector<string> > out;

	unsigned int cycle = 0;
	vector<string> split = extractBraketed(_input, _leftBraket, _rightBraket);
	while (true) {
		vector<string> tokens0 = tokenize(split[0], _delimiter);
		for (vector<string>::iterator k=tokens0.begin(); k!=tokens0.end(); k++) {
			out.push_back(vector<string>(1, *k));
		}
		vector<string> tokens1 = tokenize(split[1], _delimiter);
		if (tokens1.size() > 0) {
			out.push_back(tokens1);
		}
		split = extractBraketed(split[2], _leftBraket, _rightBraket);
		cycle++;
		if (trim(split[2]).length() == 0) {
			vector<string> tokens0 = tokenize(split[0], _delimiter);
			for (vector<string>::iterator k=tokens0.begin(); k!=tokens0.end(); k++) {
				out.push_back(vector<string>(1, *k));
			}
			vector<string> tokens1 = tokenize(split[1], _delimiter);
			if (tokens1.size() > 0) {
				out.push_back(tokens1);
			}
			break;
		}
	}
	return out;
}

vector<string> MslTools::extractBraketed(const string & _input, const string & _leftBraket, const string _rightBraket) {
	vector<string> out;
	size_t posL  = _input.find(_leftBraket);
	size_t posR  = _input.find(_rightBraket, posL);
	if (posL != std::string::npos && posR != std::string::npos) {
		out.push_back(_input.substr(0, posL));
		out.push_back(_input.substr(posL+1, posR-posL-1));
		out.push_back(_input.substr(posR+1, _input.length()-posR-1));
	} else {
		out.push_back(_input);
		out.push_back("");
		out.push_back("");
	}
	return out;
}


bool MslTools::splitIntAndString(const string & _input, int & _intResult, string & _stringResult) {

	bool foundDigit = false;

	_intResult = 0;
	_stringResult = "";

	string trimmed = trim(_input);

	int scaleFactor = 1;
	int startIndex  = 0;
	if (trimmed[0] == '-'){
		scaleFactor = -1;
		startIndex = 1;
	}
	for (int i=startIndex; i<trimmed.size(); i++) {

		int asciiCode = trimmed[i];

		if (asciiCode >= 48 && asciiCode <= 57) {
			// it is a digit
			foundDigit = true;
			int digit = trimmed[i] - '0';
			_intResult *= 10;
			_intResult += digit;
		} else {
			_stringResult = trim(trimmed.substr(i, trimmed.size() - i));
			_intResult *= scaleFactor;
			return foundDigit;
		}
	}

	_intResult *= scaleFactor;
	return foundDigit;
}

string MslTools::joinLines(const vector<string> & _input, const string & _spacer) {
	string out;
	for (vector<string>::const_iterator k=_input.begin(); k!=_input.end(); k++) {
		if (k != _input.begin()) {
			out += _spacer;
		}
		out += *k;
	}
	return out;
}

// CHECK CONTENT OF STRING
bool MslTools::isDigitChars(string _input) {
	for (int i=0; i<_input.size(); i++) {
		int asciiCode = _input[i];
		if (asciiCode < 48 && asciiCode > 57) {
			return false;
		}
	}
	return true;
}

bool MslTools::isAlphaNumericChars(string _input) {
	for (int i=0; i<_input.size(); i++) {
		int asciiCode = _input[i];
		if (!(asciiCode >= 48 && asciiCode <= 57) && !(asciiCode >= 65 && asciiCode <= 90) && !(asciiCode >= 97 && asciiCode <= 122)) {
			return false;
		}
	}
	return true;
}

bool MslTools::isAlphaChars(string _input) {
	for (int i=0; i<_input.size(); i++) {
		int asciiCode = _input[i];
		if (!(asciiCode >= 65 && asciiCode <= 90) && !(asciiCode >= 97 && asciiCode <= 122)) {
			return false;
		}
	}
	return true;
}

bool MslTools::isWhiteSpaces(string _input) {
	// space 32, tab \t 9, line feed \n 10, carriage return \r 13
	for (int i=0; i<_input.size(); i++) {
		int asciiCode = _input[i];
		if (asciiCode != 32 && asciiCode != 9 && asciiCode != 10 && asciiCode != 13) {
			return false;
		}
	}
	return true;
}

unsigned int MslTools::hamming_distance(const std::string &_str1, const std::string &_str2){
  
  if (_str1.size() != _str2.size()){
    throw MslSizeException("hamming_distance requires equal length strings");
  }

  return std::inner_product(_str1.begin(),_str1.end(), _str2.begin(), 0, std::plus<unsigned int>(), std::not2(std::equal_to<std::string::value_type>()));
  
}

vector<string> MslTools::joinConnectedLines(const vector<string> & _input, string _marker, string _spacer) {
	/******************************************************
	 *
	 *  Connect lines if they have a special terminating character
	 *  (the defauls is a backlash)
	 *
	 *  Example 
	 *
	 *    This is a first line
	 *    This is a second line
	 *    This is a third line, backslashed,\
	 *    and this is the continuation of the third line
	 *    This is a fourth line
	 *
	 *  gives:
	 *    This is a first line
	 *    This is a second line
	 *    This is a third line, backslashed, and this is the continuation of the third line
	 *    This is a fourth line
	 *******************************************************/

	vector<string> out;
	for (vector<string>::const_iterator k=_input.begin(); k!=_input.end(); k++) {
		out.push_back("");
		bool first = true;
		while (k->size() > 0 && k->substr(k->size()-1, 1) == _marker) {
			if (first) {
				first = false;
			} else {
				out.back() += _spacer;
			}
			out.back() += k->substr(0, k->size()-1);
			k++;
			if (k==_input.end()) {
				return out;
			}
		}
		if (first) {
			first = false;
		} else {
			out.back() += _spacer;
		}
		out.back() += *k;
	}
	return out;
}


string MslTools::uncomment(const string & _input, const string & _commentString) {

	size_t foundComment = _input.find(_commentString);
	if (foundComment == string::npos) {
		// no comment string, we are done
		return _input;
	}

	bool openSingleQuote = false;
	bool openDoubleQuote = false;
	for (unsigned int i=0; i< _input.size() - _commentString.size() + 1; i++) {
		// truncate a line at the # character
		if (!openSingleQuote && !openDoubleQuote && _input.substr(i, _commentString.size()) == _commentString) {
			return _input.substr(0,i);
		}
		// but do not consider if the # is within quotes
		if (!openSingleQuote && _input[i] == '"') {
			openDoubleQuote = !openDoubleQuote;
		}
		if (!openDoubleQuote && _input[i] == '\'') {
			openSingleQuote = !openSingleQuote;
		}
	}
	return _input;
}

vector<string> MslTools::uncomment(const vector<string> & _input, const string & _commentString) {
	vector<string> out;
	for (vector<string>::const_iterator k=_input.begin(); k!=_input.end(); k++) {
		out.push_back(uncomment(*k, _commentString));
	}
	return out;
}

vector<string> MslTools::removeEmptyLines(const vector<string> & _input) {
	vector<string> out;
	for (vector<string>::const_iterator k=_input.begin(); k!=_input.end(); k++) {
		// remove empty lines
		if (trim(*k).size() != 0) {
			out.push_back(*k);
		}
	}
	return out;
}



string MslTools::getFileName(string fullpath){

	string name = fullpath;
	
	// Check if fullpath only has a filename
	/*
	if (fullpath.find("/") == string::npos) {
		return name;
	}
	*/
	
	// Assume linux paths '/' , sorry windows users!
	vector<string> paths = tokenize(name, "/");

	if (paths.size() > 0){
		name = paths[paths.size()-1];
	}
	
	// If no extension, return the name
	if (name.find(".") == string::npos) {
		return name;
	}
	
	size_t pos = name.find_last_of(".");
	if (pos == std::string::npos || pos == 0) {
		return name;
	}
	
	name.erase(pos);

	return name;
}




double MslTools::correlate(vector<double> _data1, vector<double> _data2){

	// GSL Way... (GSL 1.10+ required)
	/*
	double d1[_data1.size()];
	double d2[_data1.size()];
	for (uint i = 0; i < _data1.size();i++){
		d1[i] = _data1[i];
		d2[i] = _data2[i];
	}

	// This gets the Pearson correlation coefiicient (the 1's are strides: skipping number in array)
	//  double r =  gsl_stats_correlation(d1,1,d2,1,_data1.size());  

	*/         


	// Pearson product-moment correlation coefficient 
	//                               n sum(x y) - sum(x) sum(y)
	//        r(x, y) := ---------------------------------------------------
	//                                        2                           2
	//                   sqrt(n sum(x x) - sum (x)) sqrt(n sum(y y) - sumy )


	double sumX;
	double sumY;
	double sumXY;
	double sumXX;
	double sumYY;
	sumX = sumY = sumXY = sumXX = sumYY = 0.0;
	for (uint i = 0 ; i < _data1.size();i++){

		sumX += _data1[i];
		sumY += _data2[i];
		sumXY += (_data1[i] *_data2[i]);
		sumXX += (_data1[i] *_data1[i]);
		sumYY += (_data2[i] *_data2[i]);
		
	}

	int n = _data1.size();
	double r =  (n * sumXY - sumX*sumY) / (sqrt(n * sumXX - (sumX*sumX)) * sqrt(n * sumYY - (sumY*sumY)));
	return r;
}



bool MslTools::readTextFile(vector<string> & _container, const string & _filename) {

	_container.clear();

	ifstream fs;
	fs.open(_filename.c_str());
	if (fs.fail()) {
		return false;
	}

	while (true) {
		string line;
		getline(fs,line);

		if (fs.fail()) {
			// no more lines to read from file, quit the while loop
			break;
		}

		_container.push_back(line);

	}
	fs.close();

	return true;

}

bool MslTools::getPdbsFromFile(string _filename, vector<string> &_pdbs, string _pdbPath, string _cifPath) {
	vector<string> tmp;
	return getPdbsFromFile(_filename, _pdbs, tmp, _pdbPath, _cifPath, false);
}	
bool MslTools::getPdbsFromFile(string _filename, vector<string> &_pdbs, vector<string> &_chainIds, string _pdbPath, string _cifPath, bool _parseChains) {

	ifstream fs;

	fs.open(_filename.c_str());
	if (fs.fail()){
		cerr<<"Cannot open file "<<_filename<<endl;
		exit(1);
	}

	bool isFirstValidLine = true;
	bool isPDBFile = false;

	while (true) {
		string line;
		getline(fs, line);

		if (fs.fail()) {
			// no more lines to read, quit the while.
			break;
		}
		
		// Skip if blank line or first char is a '#' 
		if (line == "" || line[0] == '#'){
			continue;
		}
	
		// Check if the line is a PDB file or a cull-PDB type file (PDBidCHAINid) for the first non-blank line
		if (isFirstValidLine) {
			if (line.size() != 8 || line.substr(line.size() - 4) != ".pdb") {
				isPDBFile = false;
				cout << "Assuming the file type is PDB chain format (PDBidCHAINid)." << endl;
			} else {
				isPDBFile = true;
				cout << "Assuming the file type is PDB file." << endl;
			}
			isFirstValidLine = false;
		}

		// Deal with non-pdb file list
		if (!isPDBFile) {

			// First check that pdbPath exists
			if (!MslTools::directoryExists(_pdbPath)) {
				cerr << "PDB path does not exist: " << _pdbPath << endl;
			}

			// Check if the file exists in the pdbPath
			string prefix = line.substr(0, 4);
			transform(prefix.begin(), prefix.end(), prefix.begin(), ::tolower);
			string filename     = _pdbPath + "/" + prefix.substr(1, 2) + "/" + prefix + ".pdb";
			string cif_filename = _cifPath + "/" + prefix.substr(1, 2)  + "/" + prefix + ".cif";

			// If file exists, add to pdbs and chain ids
			if (MslTools::fileExists(filename)) {
				
				_pdbs.push_back(filename);

				if (_parseChains) {
					// get chainid from line character 5 until a whitespace
					string chainId = line.substr(4, line.find(" ", 4) - 4);
					_chainIds.push_back(chainId);
				}

			} else {

				// Try to check if the file exists in the cifPath
				bool tryRemote = true;
				if (MslTools::directoryExists(_cifPath)) {

					if (MslTools::fileExists(cif_filename)){
						_pdbs.push_back(cif_filename);
						if (_parseChains) {
							string chainId = line.substr(4, line.find(" ", 4) - 4);
							_chainIds.push_back(chainId);
						}
		
						tryRemote = false;
					} else {
						cerr << "Cannot find PDB file: " << filename << " or CIF file: " << cif_filename << " trying remote." << endl;
					}

				} else {
					cerr << "CIF path does not exist: " << _cifPath << " trying remote." << endl;
				}
				
					
				if (tryRemote) {
					// Get from remote server
					stringstream ss;
					ss << MslTools::toUpper(prefix) << ".cif";
					string filename2 = ss.str();

					string pdbUrl = "https://files.rcsb.org/download/" + filename2;
					
					string command = "wget -O " + filename2 + " " + pdbUrl;
					int retVal = system(command.c_str());
					
					if (MslTools::fileExists(filename2)){
						_pdbs.push_back(filename2);

						if (_parseChains) {
							string chainId = line.substr(4, line.find(" ", 4) - 4);
							_chainIds.push_back(chainId);
						}
					} else {
						cerr << "SKIPPING FILE: Cannot find get remote get file: "<<filename2<<endl;
					}
				}
			}	
		} else {
			_pdbs.push_back(line);
		}
	}

	return true;

}

bool MslTools::fileExists(string _filename){
  // This way is broken GCC 9.2.0
  //ifstream afile(_filename.c_str());
  //return afile.fail();
  return (access(_filename.c_str(), F_OK) != -1);
}

// a function that tests if a directory exists
bool MslTools::directoryExists(std::string _dir) {
	struct stat info;
	if (stat(_dir.c_str(), &info) != 0) {
		return false;
	}
	return (info.st_mode & S_IFDIR);
}

string MslTools::pathRoot(string _path) {
	// /the/path/theFile.ext -> /the/path/theFile
	for (int i=_path.size()-1; i>=0; i--) {
		if (_path.substr(i,1) == "/") {
			return _path;
		} else if (_path.substr(i,1) == ".") {
			return _path.substr(0,i);
		}
	}
	return _path;
}

string MslTools::pathHead(string _path) {
	// /the/path/theFile.ext -> /the/path
	for (int i=_path.size()-1; i>=0; i--) {
		if (_path.substr(i,1) == "/") {
			return _path.substr(0,i);
		}
	}
	return _path;
}

string MslTools::pathExtension(string _path) {
	// /the/path/theFile.ext -> ext
	for (int i=_path.size()-1; i>=0; i--) {
		if (_path.substr(i,1) == "/") {
			return "";
		} else if (_path.substr(i,1) == ".") {
			return _path.substr(i+1,_path.size()-i-1);
		}
	}
	return "";
}

string MslTools::pathTail(string _path) {
	// /the/path/theFile.ext -> theFile.ext
	for (int i=_path.size()-1; i>=0; i--) {
		if (_path.substr(i,1) == "/") {
			return _path.substr(i+1,_path.size()-i-1);
		}
	}
	return _path;
}




double MslTools::smartRound(double coord, double gridSize){

  if (gridSize > 1) {
    int cint = int(coord*100);
    int gint = int(gridSize *100);
    double remainder = gint - (cint % gint);
    if (remainder / (double)gint < 0.5)
      return int(coord + remainder/100);

    return int((coord*100 - (gridSize*100 - remainder))/100);

  } else {

    double frac_coord       = coord - int(coord);
    int shift_frac_coord = int(100 * frac_coord);
    double correction = int(shift_frac_coord % int(gridSize *100));
    double frac_result = (shift_frac_coord - correction)/100;
    if (correction/(int (gridSize*100)) >= 0.5){
      frac_result += gridSize;
      if (frac_result > 1) frac_result = int(frac_result);	
    }

    return int(coord) + (frac_result);
  }
}


double MslTools::mod(double x, double y){
	return x-y *floor(x/y);
}


Real MslTools::round(Real value) {
    Real half = (Real) 0.5f;
    int val = (int)( (value > (Real)0.0f)?(value + half):(value - half) );
    return( (Real)val );
}


string MslTools::createDir(string _name){

	// First make sure name is unique
	string prevOutputdir = _name;
	string outputdir = outputFileNameParser(_name);
	while (true) {
		if (prevOutputdir == outputdir) {
			// the outputdir doesn't have randomly generated parts
			break;
		}
		ifstream test_fs;
		test_fs.open(outputdir.c_str());
		if (test_fs.is_open()) {
			// directory exists, try again
			prevOutputdir = outputdir;
			outputdir = outputFileNameParser(outputdir);
			test_fs.close();
		} else {
			test_fs.close();
			break;
		}
	}

	// Now make the directory
	mkNestedDir(outputdir, 0755); 

	return outputdir;
}

string MslTools::outputFileNameParser(string _name) {
	vector<unsigned int> startSub;
	vector<unsigned int> endSub;
	vector<unsigned int> size;
	vector<string> type;
	for (unsigned int i=0; i<_name.size(); i++) {
		if (_name.substr(i, 1) == "%") {
			unsigned int start = i;
			i++;
			bool invalid = false;
			string digits;
			while (i<_name.size() && _name.substr(i,1) != "s" && _name.substr(i,1) != "r" && _name.substr(i,1) != "t" && _name.substr(i,1) != "d") {
				if (_name[i] > 47 && _name[i] < 58) {
					digits += _name.substr(i,1);
					i++;
				} else {
					invalid = true;
					break;
				}
			}
			if (!invalid) {
				if (digits.size() == 0) {
					digits = "1";
				}
				unsigned int number = toInt(digits);
				if ( _name.substr(i,1) != "s" || _name.substr(i,1) != "r" || _name.substr(i,1) != "t" || _name.substr(i,1) != "d") {
					startSub.push_back(start);
					endSub.push_back(i);
					size.push_back(number);
					type.push_back(_name.substr(i,1));
				}
			}
		}
	}

	time_t currentTime = time(NULL);; 
	tm *curr = localtime(&currentTime);
	//time(&currentTime);
	stringstream ss;
	ss << (int)currentTime;
	string timeString = ss.str();
	char c [1000];
	sprintf(c, "%04d%02d%02d", curr->tm_year+1900, curr->tm_mon+1, curr->tm_mday);
	string dateStamp = c;
	sprintf(c, "%02d%02d%02d", curr->tm_hour, curr->tm_min, curr->tm_sec);
	string timeStamp = c;

	for (int i=startSub.size()-1; i>=0; i--) {
		if (type[i] == "s") {
			while (timeString.size() < size[i]) {
				timeString = "0" + timeString;
			}
			if (timeString.size() > size[i]) {
				timeString = timeString.substr(timeString.size() - size[i], size[i]);
			}
			_name.replace(startSub[i], endSub[i] - startSub[i] + 1, timeString);
		} else if (type[i] == "d") {
			_name.replace(startSub[i], endSub[i] - startSub[i] + 1, dateStamp);
		} else if (type[i] == "t") {
			_name.replace(startSub[i], endSub[i] - startSub[i] + 1, timeStamp);
		} else if (type[i] == "r") {
			_name.replace(startSub[i], endSub[i] - startSub[i] + 1, getRandomAlphaNumString(size[i],true));
		} 
	}
	return _name;

}

bool MslTools::mkNestedDir(string _dir, mode_t _mode) {

	_dir = trim(_dir);
	vector<string> directories;
	vector<string>::iterator k;

	string prevPath = "";
	if (_dir.substr(0,1) == "/") {
		prevPath = "/";
		_dir = _dir.substr(1,_dir.size()-1);
	}
	
	bool open = false;
	for (int i=0; i< _dir.size(); i++) {
		if (_dir.substr(i,1) == "/") {
			open = false;
			continue;
		} else {
			if (open) {
				*k += _dir[i];
			} else {
				open = true;
				directories.push_back(_dir.substr(i,1));
				k = directories.end() - 1;
			}
		}
	}

	for (unsigned int i=0; i<directories.size(); i++) {
		// don't make . and ..
		if (directories[i] != "." || directories[i] != "..") {
			string dirToMake = prevPath + directories[i];
			
			if (mkdir(dirToMake.c_str(), _mode) != 0) {
				// failed because it did already exist?
				if (errno != EEXIST) {
					// no, we failed
					return false;
				}
			}
			prevPath += directories[i] + (string)"/";
		}
	}
	return true;
}

string MslTools::getRandomAlphaNumString(unsigned int _size, bool _alphaOnly) {
	vector<string> characters;
	characters.push_back("a");
	characters.push_back("b");
	characters.push_back("c");
	characters.push_back("d");
	characters.push_back("e");
	characters.push_back("f");
	characters.push_back("g");
	characters.push_back("h");
	characters.push_back("i");
	characters.push_back("j");
	characters.push_back("k");
	characters.push_back("l");
	characters.push_back("m");
	characters.push_back("n");
	characters.push_back("o");
	characters.push_back("p");
	characters.push_back("q");
	characters.push_back("r");
	characters.push_back("s");
	characters.push_back("t");
	characters.push_back("u");
	characters.push_back("v");
	characters.push_back("w");
	characters.push_back("x");
	characters.push_back("y");
	characters.push_back("z");
	characters.push_back("A");
	characters.push_back("B");
	characters.push_back("C");
	characters.push_back("D");
	characters.push_back("E");
	characters.push_back("F");
	characters.push_back("G");
	characters.push_back("H");
	characters.push_back("I");
	characters.push_back("J");
	characters.push_back("K");
	characters.push_back("L");
	characters.push_back("M");
	characters.push_back("N");
	characters.push_back("O");
	characters.push_back("P");
	characters.push_back("Q");
	characters.push_back("R");
	characters.push_back("S");
	characters.push_back("T");
	characters.push_back("U");
	characters.push_back("V");
	characters.push_back("W");
	characters.push_back("X");
	characters.push_back("Y");
	characters.push_back("Z");
	characters.push_back("0");
	characters.push_back("1");
	characters.push_back("2");
	characters.push_back("3");
	characters.push_back("4");
	characters.push_back("5");
	characters.push_back("6");
	characters.push_back("7");
	characters.push_back("8");
	characters.push_back("9");
	string out;

	int size = characters.size();
	if (_alphaOnly){
		size = 52; //  53 ?
	}
/*

#ifdef __GSL__
	RandomNumberGenerator rng;
	rng.setRNGType("knuthran2");
	rng.setRNGTimeBasedSeed();
#endif

	for (unsigned int i=0; i<_size; i++) {
#ifdef __GSL__
		int randomN = (int)(rng.getRandomInt() % size);
#else
		//int randomN = (int)(-(double)rand()/(double)(RAND_MAX+1) * (double)characters.size());
		int randomN = rand() % size;
#endif
		out += characters[randomN];
	}
*/	
	
	RandomNumberGenerator rng;

	for (unsigned int i=0; i<_size; i++) {
		int randomN = (int)(rng.getRandomInt(size-1));
		out += characters[randomN];
	}
	
	return out;

}

/*
unsigned int MslTools::getRandomInt(unsigned int _max) {
#ifdef __GSL__
	RandomNumberGenerator rng;
	rng.setRNGType("knuthran2");
	rng.setRNGTimeBasedSeed();
	return (int)(rng.getRandomInt() % _max);
#else
	return rand() % _max;
#endif
}
*/

void MslTools::loadAAConversionTables() {
	threeToOneLetter["ALA"] = "A";
	oneToThreeLetter["A"]   = "ALA";

	threeToOneLetter["CYS"] = "C";
	oneToThreeLetter["C"]   = "CYS";

	threeToOneLetter["ASP"] = "D";
	oneToThreeLetter["D"]   = "ASP";

	threeToOneLetter["GLU"] = "E";
	oneToThreeLetter["E"]   = "GLU";

	threeToOneLetter["PHE"] = "F";
	oneToThreeLetter["F"]   = "PHE";

	threeToOneLetter["GLY"] = "G";
	oneToThreeLetter["G"]   = "GLY";

	threeToOneLetter["HIS"] = "H";
	oneToThreeLetter["H"]   = "HIS";

	threeToOneLetter["HSE"] = "H";
	threeToOneLetter["HSP"] = "H";
	threeToOneLetter["HSD"] = "H";

	threeToOneLetter["ILE"] = "I";
	oneToThreeLetter["I"]   = "ILE";

	threeToOneLetter["LYS"] = "K";
	oneToThreeLetter["K"]   = "LYS";

	threeToOneLetter["LEU"] = "L";
	oneToThreeLetter["L"]   = "LEU";

	threeToOneLetter["MET"] = "M";
	threeToOneLetter["MSE"] = "M";
	oneToThreeLetter["M"]   = "MET";

	threeToOneLetter["ASN"] = "N";
	oneToThreeLetter["N"]   = "ASN";

	threeToOneLetter["PRO"] = "P";
	oneToThreeLetter["P"]   = "PRO";

	threeToOneLetter["GLN"] = "Q";
	oneToThreeLetter["Q"]   = "GLN";

	threeToOneLetter["ARG"] = "R";
	oneToThreeLetter["R"]   = "ARG";

	threeToOneLetter["SER"] = "S";
	oneToThreeLetter["S"]   = "SER";

	threeToOneLetter["THR"] = "T";
	oneToThreeLetter["T"]   = "THR";

	threeToOneLetter["VAL"] = "V";
	oneToThreeLetter["V"]   = "VAL";

	threeToOneLetter["TRP"] = "W";
	oneToThreeLetter["W"]   = "TRP";

	threeToOneLetter["TYR"] = "Y";
	oneToThreeLetter["Y"]   = "TYR";
}
string MslTools::getOneLetterCode(string threeLetterCode) { 
    string oneLetterCode = "X";
    if (threeToOneLetter.size()<= 0) 
        loadAAConversionTables();
    threeLetterCode = toUpper(threeLetterCode);
    map<string, string>::iterator entry = threeToOneLetter.find(threeLetterCode);
    if (entry != threeToOneLetter.end())
        oneLetterCode = entry->second;
    return oneLetterCode;
}

string MslTools::getThreeLetterCode(string oneLetterCode) { 
    string threeLetterCode = "XXX";
    if (oneToThreeLetter.size()<= 0) 
        loadAAConversionTables();
    oneLetterCode = toUpper(oneLetterCode);
    map<string, string>::iterator entry = oneToThreeLetter.find(oneLetterCode);
    if (entry != oneToThreeLetter.end())
        threeLetterCode = entry->second;
    return threeLetterCode;
}


bool MslTools::sortPairIntDoubleAscending(const pair<int,double> &left, const pair<int,double> &right){
	return left.second < right.second;
}
bool MslTools::sortPairIntDoubleDecending(const pair<int,double> &left, const pair<int,double> &right){
	return left.second > right.second;
}



double MslTools::setPrecision(double _d, unsigned int _significantDigits) {
	return (double)int(_d * pow(10.0, (int)_significantDigits) + 0.5) / pow(10.0, (int)_significantDigits);
}




bool MslTools::sortByResnumIcodeAscending(int _resnum1, string _icode1, int _resnum2, string _icode2) {
	if (_resnum1 < _resnum2 || (_resnum1 == _resnum2 && _icode1.c_str()[0] < _icode1.c_str()[1])) {
		return true;
	}
	return false;
}

double MslTools::partition(vector<double> & _vec, unsigned int _left, unsigned int _right, unsigned int _pivotIndex) {
	// used by quickSort
	double tmpPivotValue = _vec[_pivotIndex];
	_vec[_pivotIndex] = _vec[_right];
	_vec[_right] = tmpPivotValue;
	unsigned int storeIndex = _left;
	double tmp = 0.0;
	for (unsigned int i=_left; i< _right; i++) {
		if (_vec[i] <= tmpPivotValue) {
			if (i != storeIndex) {
				tmp = _vec[i];
				_vec[i] = _vec[storeIndex];
				_vec[storeIndex] = tmp;
			}
			storeIndex++;
		}
	}
	tmp = _vec[_right];
	_vec[_right] = _vec[storeIndex];
	_vec[storeIndex] = tmp;
	return storeIndex;
}

void MslTools::quickSort(vector<double> & _vec, unsigned int _left, unsigned int _right) {
	if (_right-_left < 1) {
		return;
	}
	unsigned int pivotIndex = int(partition(_vec, _left, _right, _left));
	if (pivotIndex > _left) {
		quickSort(_vec, _left, pivotIndex-1);
	}
	if (pivotIndex < _right) {
		quickSort(_vec, pivotIndex+1, _right);
	}
}

void MslTools::quickSort(vector<double> & _vec) {
	if (_vec.size() < 2) {
		return;
	}
	quickSort(_vec, 0, _vec.size()-1);
}

double MslTools::partitionWithIndex(vector<double> & _vec, unsigned int _left, unsigned int _right, unsigned int _pivotIndex, vector<unsigned int> &_index) {
	// used by quickSortWithIndex
	double tmpPivotValue = _vec[_pivotIndex];
	_vec[_pivotIndex] = _vec[_right];
	_vec[_right] = tmpPivotValue;
	unsigned int tmpPivotIndex = _index[_pivotIndex];
	_index[_pivotIndex] = _index[_right];
	_index[_right] = tmpPivotIndex;
	unsigned int storeIndex = _left;
	unsigned int tmpI = 0;
	double tmp = 0;
	for (unsigned int i=_left; i< _right; i++) {
		if (_vec[i] <= tmpPivotValue) {
			if (i != storeIndex) {
				tmp = _vec[i];
				_vec[i] = _vec[storeIndex];
				_vec[storeIndex] = tmp;
				tmpI = _index[i];
				_index[i] = _index[storeIndex];
				_index[storeIndex] = tmpI;
			}
			storeIndex++;
		}
	}
	tmp = _vec[_right];
	_vec[_right] = _vec[storeIndex];
	_vec[storeIndex] = tmp;
	tmpI = _index[_right];
	_index[_right] = _index[storeIndex];
	_index[storeIndex] = tmpI;
	return storeIndex;
}

void MslTools::quickSortWithIndex(vector<double> & _vec, unsigned int _left, unsigned int _right, vector<unsigned int> &_index) {
	if (_vec.size() != _index.size()) {
		cerr << "ERROR 31286: number of element mismatch _vec: " << _vec.size() << " != _index: " << _index.size() << " at void MslTools::quickSortWithIndex(vector<double> & _vec, unsigned int _left, unsigned int _right, vector<unsigned int> &_index)" << endl;
		exit(31286);
	}
	if (_right-_left < 1) {
		return;
	}
	unsigned int pivotIndex = int(partitionWithIndex(_vec, _left, _right, _left, _index));
	if (pivotIndex > _left) {
		quickSortWithIndex(_vec, _left, pivotIndex-1, _index);
	}
	if (pivotIndex < _right) {
		quickSortWithIndex(_vec, pivotIndex+1, _right, _index);
	}
}

void MslTools::quickSortWithIndex(vector<double> & _vec, vector<unsigned int> &_index) {
	if (_vec.size() < 2) {
		return;
	}
	quickSortWithIndex(_vec, 0, _vec.size()-1, _index);
}

void MslTools::normalizeVector(vector<double> & _vec) {
	// normalize a vector so that the sum is = 1
	if (_vec.size() == 0) {
		cerr << "WARNING 4843: vector size is zero in void cubTools::normalizeCumulativeVector(vector<double> & _vec)" << endl;
		return;
	}
	double sum = 0.0;
	for (unsigned int i=0; i<_vec.size(); i++) {
		sum += _vec[i];
	}
	if (sum == 0.0) {
		cerr << "ERROR 4846: cannot renormalize by zero total in void cubTools::normalizeVector(vector<double> & _vec)" << endl;
		exit(4846);
	}
	for (unsigned int i=0; i<_vec.size(); i++) {
		_vec[i] /= sum;
	}

}

void MslTools::normalizeCumulativeVector(vector<double> & _vec) {
	// normalize a cumulative vector so that the last element is = 1.
	// since it element is supposed to be the sum of the previous
	// the function issues warnings if an element is smaller than the 
	// previous
	if (_vec.size() == 0) {
		cerr << "ERROR 4850: vector size is zero in void cubTools::normalizeCumulativeVector(vector<double> & _vec)" << endl;
		exit(4850);
	} else {
		if (_vec.back() == 0.0) {
			cerr << "ERROR 4853: cannot renormalize by zeroed last element in void cubTools::normalizeCumulativeVector(vector<double> & _vec)" << endl;
			exit(4853);
		}
	}
	double prev = 0.0;
	for (unsigned int i=0; i<_vec.size(); i++) {
		if (_vec[i] < prev) {
			cerr << "WARNING 4854: element " << i << " (" << _vec[i] << " is less than previous element (" << prev << " in void cubTools::normalizeCumulativeVector(vector<double> & _vec)" << endl;
		}
		prev = _vec[i];
		_vec[i] /= _vec.back();
	}
}


void MslTools::rgb2hsv(vector<double> &_rgb, vector<double> &_hsv){
	_hsv.clear();
	_hsv.push_back(0.0);
	_hsv.push_back(0.0);
	_hsv.push_back(0.0);

	if (_rgb.size() != 3){
		cerr << "ERROR MslTools::rgb2hsv() _rgb not of size 3\n";
		exit(1);
	}


	double max = *(std::max_element(_rgb.begin(),_rgb.end()));
	double min = *(std::min_element(_rgb.begin(),_rgb.end()));
	
	_hsv[2] = max;
	_hsv[1] = (max==0)?0.0f:(1.0f - min/max);
	if(max == min)
		_hsv[0] = 0;
	else if(max == _rgb[0]) {
	        float tmp = (_rgb[1]-_rgb[2])/(max-min);
		_hsv[0] = fmod(60.0f *  tmp + 360.0f, 360.0f);
	} else if(max == _rgb[1])
		_hsv[0] = (60.0f * (_rgb[2]-_rgb[0])/(max-min) + 120.0f);
	else
		_hsv[0] = (60.0f * (_rgb[0]-_rgb[1])/(max-min) + 240.0f);
    
}
void MslTools::hsv2rgb(vector<double> &_hsv, vector<double> &_rgb){

	_rgb.clear();
	_rgb.push_back(0.0);
	_rgb.push_back(0.0);
	_rgb.push_back(0.0);
	if (_hsv.size() != 3){
		cerr << "ERROR MslTools::hsv2rgb() _hsv not of size 3\n";
		exit(1);
	}

    int hi = (int)(_hsv[0]/60.0f);
    float f = _hsv[0]/60.0f - hi;
    hi = int(mod(double(hi), 6.0f));
    
    float p = _hsv[2] * (1.0f - _hsv[1]);
    float q = _hsv[2] * (1.0f - f*_hsv[1]);
    float t = _hsv[2] * (1.0f - (1.0f-f)*_hsv[1]);
    switch(hi) {
        case 0:
            _rgb[0] = _hsv[2]; _rgb[1] = t; _rgb[2] = p;
            break;
        case 1:
            _rgb[0] = q; _rgb[1] = _hsv[2]; _rgb[2] = p;
            break;
        case 2:
            _rgb[0] = p; _rgb[1] = _hsv[2]; _rgb[2] = t;
            break;
        case 3:
            _rgb[0] = p; _rgb[1] = q; _rgb[2] = _hsv[2];
            break;
        case 4:
            _rgb[0] = t; _rgb[1] = p; _rgb[2] = _hsv[2];
            break;
        case 5:
            _rgb[0] = _hsv[2]; _rgb[1] = p; _rgb[2] = q;
            break;
        default:
	  cerr << "ERROR1 MslTools::hsv2rgb() problem. see code("<<hi<<")\n";
		exit(1);
		break;
    }
    
}

vector<double> MslTools::getRGB(vector<double> &_startRGB, vector <double> &_endRGB, double _minValue, double _maxValue, double _value){
	
	// Bound value to one of our limits
	if (_value > _maxValue) _value = _maxValue;
	if (_value < _minValue) _value = _minValue;

	// Convert to HSV
	vector<double> hsvStart(3,0.0);
	rgb2hsv(_startRGB,hsvStart);

	vector<double> hsvEnd(3,0.0);
	rgb2hsv(_endRGB,hsvEnd);

	// Convert value from 0 to 1.
	double normVal = (_value - _minValue) / (_maxValue - _minValue);

	
	// Linear interpolation between start and end for each H,S,V values
	vector<double> resultHSV(3,0.0);
	for (uint i = 1; i < 3; i++){
		resultHSV[i] = (hsvStart[i] * normVal) + (hsvEnd[i] * 1-normVal);
	}
	
	float tmp = (hsvEnd[0] - hsvStart[0]);
	double hdiff =  fmod(tmp, 360.f);
	if(hdiff > 180.0f) {
		hdiff = 360.0f-hdiff;
		normVal = -normVal;
	}
	tmp = hsvStart[0] + hdiff * normVal;
	resultHSV[0] = fmod(tmp, 360.0f);

	vector<double> resultRGB(3,0.0);
	hsv2rgb(resultHSV,resultRGB);
	return resultRGB;
}


string MslTools::getMSLversion() {
	return (string) "MSL v." + (string)MSLVERSION + (string)" of " + (string)MSLDATE;
}

// The Atom Id is in the form of "A 37 CA" or "A 37A CA" with an insertion code
string MslTools::getAtomId(string _chainid, int _resnum, string _icode, string _atomName, unsigned int _skiplevel) {
	_chainid = MslTools::trim(_chainid);
	_icode = MslTools::trim(_icode);
	_atomName = MslTools::trim(_atomName);
	char c [1000];
	if (_skiplevel == 0) {
		sprintf(c, "%s,%d%s,%s", _chainid.c_str(), _resnum, _icode.c_str(), _atomName.c_str());
	} else if (_skiplevel == 1) {
		sprintf(c, "%d%s,%s", _resnum, _icode.c_str(), _atomName.c_str());
	} else {
		sprintf(c, "%s", _atomName.c_str());
	}
	return (string)c;
}
/*
string MslTools::getAtomId(string _chainid, int _resnum, string _atomName) {
	return getAtomId(_chainid, _resnum, "", _atomName);
}
*/
bool MslTools::parseAtomId(string _atomId, string & _chainid, int & _resnum, string & _icode, string & _atomName, unsigned int _skiplevels) {
	// parses "A 37 CA" or "A 37A CA" with icode
	// a skip level of 1 allows to pass "37 CA" (for calls from chain)
	// a skip level of 2 allows to pass "CA" (for calls from residue or position)

	vector<string> tokens = MslTools::tokenize(_atomId, ",", true);
	if (tokens.size() == 1) {
		// no comma in string
		tokens = MslTools::tokenize( _atomId, " ");
		if (tokens.size() == 1) {
			// no space in string
			tokens = MslTools::tokenize( _atomId, "_");
		}
	}
	_chainid = "";
	_resnum = 0;
	_icode = "";
	_atomName = "";
	if (tokens.size() != 3) {
		// not "A 37 CA" format
		if (tokens.size() + _skiplevels < 3) {
			//OK = false;
			return false;
		} else {
			if (tokens.size() == 2 && _skiplevels >= 1) {
				// "37 CA": add a blank chain
				tokens.insert(tokens.begin(), "");
			} else if (tokens.size() == 1 && _skiplevels == 2) {
				// "CA": add a blank chain a residue number
				tokens.insert(tokens.begin(), "0");
				tokens.insert(tokens.begin(), "");
			} else {
				return false;
			}
		}
	}
	for (vector<string>::iterator k=tokens.begin(); k!=tokens.end(); k++) {
		*k = MslTools::trim(*k);
	}
	_chainid = tokens[0];
	bool OK = MslTools::splitIntAndString(tokens[1], _resnum, _icode);
	_atomName = tokens[2];

	return OK;
	//_OK = true;
}

bool MslTools::compareAtomIds(string _id1, string _id2, unsigned int _skiplevels) {
	// check if two atoms ids are the same (even if they are separated by different
	// characters, it is not just a string compare, i.e. "A,37,CA" vs "A 37 CA"
	string chainid1;
	int resnum1;
	string icode1;
	string atomName1;
	string chainid2;
	int resnum2;
	string icode2;
	string atomName2;
	if (!parseAtomId(_id1, chainid1, resnum1, icode1, atomName1, _skiplevels)) {
		return false;
	}
	if (!parseAtomId(_id2, chainid2, resnum2, icode2, atomName2, _skiplevels)) {
		return false;
	}
	if (chainid1 != chainid2 || resnum1 != resnum2 || icode1 != icode2 || atomName1 != atomName2) {
		return false;
	}
	return true;

}

// The Residue Id is in the form of "A 37" or "A 37A" with an insertion code
string MslTools::getPositionId(string _chainid, int _resnum, string _icode, unsigned int _skiplevel) {
	_chainid = MslTools::trim(_chainid);
	_icode = MslTools::trim(_icode);
	char c [1000];
	if (_skiplevel == 0) {
		sprintf(c, "%s,%d%s", _chainid.c_str(), _resnum, _icode.c_str());
	} else {
		sprintf(c, "%d%s", _resnum, _icode.c_str());
	}
	return (string)c;
}
/*
string MslTools::getPositionId(string _chainid, int _resnum) {
	return getPositionId(_chainid, _resnum, "");
}
*/
bool MslTools::parsePositionId(string _posId, string & _chainid, int & _resnum, string & _icode, unsigned int _skiplevels) {
	// parses "A 37" or "A 37A" with icode
	// a skip level of 1 allows to pass "37" (for calls from chain)

	vector<string> tokens = MslTools::tokenize(_posId, ",", true);
	if (tokens.size() == 1) {
		// no comma
		tokens = MslTools::tokenize( _posId, " ");
		if (tokens.size() == 1) {
			// no comma
			tokens = MslTools::tokenize( _posId, "_");
		}
	}
	_chainid = "";
	_resnum = 0;
	_icode = "";
	if (tokens.size() != 2) {
		// not "A 37" format
		if (tokens.size() + _skiplevels < 2) {
			return false;
		} else {
			if (tokens.size() == 1 && _skiplevels == 1) {
				// "37": add a blank chain
				tokens.insert(tokens.begin(), "");
			} else {
				return false;
			}
		}
	}
	for (vector<string>::iterator k=tokens.begin(); k!=tokens.end(); k++) {
		*k = MslTools::trim(*k);
	}
	_chainid = tokens[0];
	bool OK = MslTools::splitIntAndString(tokens[1], _resnum, _icode);
	//_OK = true;
	return OK;
}
bool MslTools::compareIdentityIds(string _id1, string _id2, unsigned int _skiplevels) {
	// check if two atoms ids are the same (even if they are separated by different
	// characters, it is not just a string compare, i.e. "A,37,ILE" vs "A 37 ILE"
	string chainid1;
	int resnum1;
	string icode1;
	string identity1;
	string chainid2;
	int resnum2;
	string icode2;
	string identity2;
	if (!parseIdentityId(_id1, chainid1, resnum1, icode1, identity1, _skiplevels)) {
		return false;
	}
	if (!parseIdentityId(_id2, chainid2, resnum2, icode2, identity2, _skiplevels)) {
		return false;
	}
	if (chainid1 != chainid2 || resnum1 != resnum2 || icode1 != icode2 || identity1 != identity2) {
		return false;
	}
	return true;
}

bool MslTools::comparePositionIds(string _id1, string _id2, unsigned int _skiplevels) {
	// check if two position ids are the same (even if they are separated by different
	// characters, it is not just a string compare, i.e. "A,37" vs "A 37"
	string chainid1;
	int resnum1;
	string icode1;
	string chainid2;
	int resnum2;
	string icode2;
	if (!parsePositionId(_id1, chainid1, resnum1, icode1, _skiplevels)) {
		return false;
	}
	if (!parsePositionId(_id2, chainid2, resnum2, icode2, _skiplevels)) {
		return false;
	}
	if (chainid1 != chainid2 || resnum1 != resnum2 || icode1 != icode2) {
		return false;
	}
	return true;
}


// The Identity Id is in the form of "A 37 ILE" or "A 37A ILE" with an insertion code
string MslTools::getIdentityId(string _chainid, int _resnum, string _icode, string _identity, unsigned int _skiplevel) {
//	_chainid = MslTools::trim(_chainid);
//	_icode = MslTools::trim(_icode);
//	_identity = MslTools::trim(_identity);
//	char c [1000];
//	sprintf(c, "%s %d%s %s", _chainid.c_str(), _resnum, _icode.c_str(), _identity.c_str());
//	return (string)c;
	return getAtomId(_chainid, _resnum, _icode, _identity, _skiplevel);
}
/*
string MslTools::getIdentityId(string _chainid, int _resnum, string _identity) {
	return getAtomId(_chainid, _resnum, "", _identity);
}
*/
bool MslTools::parseIdentityId(string _residueId, string & _chainid, int & _resnum, string & _icode, string & _identity, unsigned int _skiplevels) {
	/*
	vector<string> tokens = MslTools::tokenize(_residueId, ",", true);
	if (tokens.size() != 3) {
		tokens = MslTools::tokenize( _residueId, " ");
	}
	if (tokens.size() != 3) {
		_chainid = "";
		_resnum = 0;
		_icode = "";
		_identity = "";
		_OK = false;
		return;
	}
	for (vector<string>::iterator k=tokens.begin(); k!=tokens.end(); k++) {
		*k = MslTools::trim(*k);
	}
	_chainid = tokens[0];
	MslTools::splitIntAndString(tokens[1], _resnum, _icode);
	_identity = tokens[2];
	_OK = true;
	*/
	return parseAtomId(_residueId, _chainid, _resnum, _icode, _identity, _skiplevels);
}

bool MslTools::parseMutationId(string _mutationId, string & _chainid, int & _resnum, string & _icode, string & _identity, string &_newIdentity, unsigned int _skiplevels) {

	// parses "A 37 LEU ARG" or "A 37A LEU ARG" with icode
        vector<string> tokens = MslTools::tokenize(_mutationId, ",");
	if (tokens.size() == 1) {
		// no comma in string
		tokens = MslTools::tokenize( _mutationId, " ");
		if (tokens.size() == 1) {
			// no space in string
			tokens = MslTools::tokenize( _mutationId, "_");
			
			if (tokens.size() == 1){
				cerr << "Can not tokenize string in parseMutationId("<<_mutationId<<")\n";
			}
		}
	}

	if (tokens.size() != 4) {
		return false;
	}

	for (vector<string>::iterator k=tokens.begin(); k!=tokens.end(); k++) {
		*k = MslTools::trim(*k);
	}

	_chainid = tokens[0];
	bool OK = MslTools::splitIntAndString(tokens[1], _resnum, _icode);
	if (!OK) {
		return false;
	}

	_identity = tokens[2];
	_newIdentity = tokens[3];
	return true;
}


// The Rotamer Id is in the form of "A 37 ILE 3" or "A 37A ILE 3" with an insertion code
std::string MslTools::getRotamerId(std::string _chainid, int _resnum, std::string _icode, std::string _identity, unsigned int _conformation){

	stringstream ss;
	ss << _chainid << " "<<_resnum<<_icode<<" "<<_identity<<" "<<_conformation;
	return (ss.str());
}

bool MslTools::parseRotamerId(std::string _rotamerId, std::string & _chainid, int & _resnum, std::string & _icode, std::string & _identity, unsigned int &_conformation){

	// parses "A 37 LEU 4" or "A 37A LEU 4" with icode
	vector<string> tokens = MslTools::tokenize(_rotamerId);
	if (tokens.size() == 1) {
		// no comma in string
		tokens = MslTools::tokenize( _rotamerId, " ");
		if (tokens.size() == 1) {
			// no space in string
			tokens = MslTools::tokenize( _rotamerId, "_");
			
			if (tokens.size() == 1){
				cerr << "Can not tokenize string in parseRotamerId("<<_rotamerId<<")\n";
			}
		}
	}

	if (tokens.size() != 4) {
		return false;
	}

	for (vector<string>::iterator k=tokens.begin(); k!=tokens.end(); k++) {
		*k = MslTools::trim(*k);
	}

	_chainid = tokens[0];
	bool OK = MslTools::splitIntAndString(tokens[1], _resnum, _icode);
	if (!OK) {
		return false;
	}

	_identity = tokens[2];
	_conformation = (unsigned int)MslTools::toInt(tokens[3], "Error parsing rotamer conformation into an int value from MslTools::parseRotamerId()");

	return true;

}


// The Atom of Identity Id is in the form of "A 37 ILE CA" or "A 37A ILE CA" with an insertion code
string MslTools::getAtomOfIdentityId(string _chainid, int _resnum, string _icode, string _identity, string _atomName, unsigned int _skiplevel) {
	_chainid = MslTools::trim(_chainid);
	_icode = MslTools::trim(_icode);
	_identity = MslTools::trim(_identity);
	_atomName = MslTools::trim(_atomName);
	char c [1000];
	if (_skiplevel == 0) {
		sprintf(c, "%s,%d%s,%s,%s", _chainid.c_str(), _resnum, _icode.c_str(), _identity.c_str(), _atomName.c_str());
	} else if (_skiplevel == 1) {
		sprintf(c, "%d%s,%s,%s", _resnum, _icode.c_str(), _identity.c_str(), _atomName.c_str());
	} else if (_skiplevel == 2) {
		sprintf(c, "%s,%s", _identity.c_str(), _atomName.c_str());
	} else {
		sprintf(c, "%s", _atomName.c_str());
	}
	return (string)c;
}
/*
string MslTools::getAtomOfIdentityId(string _chainid, int _resnum, string _identity, string _atomName) {
	return getAtomOfIdentityId(_chainid, _resnum, "", _identity, _atomName);
}
*/
bool MslTools::parseAtomOfIdentityId(string _atomId, string & _chainid, int & _resnum, string & _icode, string & _identity, string & _atomName, unsigned int _skiplevels) {
	// parses "A 37 ILE CA" or "A 37A ILE CA" with icode
	// a skip level of 1 allows to pass "37 ILE CA" (for calls from chain)
	// a skip level of 2 allows to pass "ILE CA" (for calls from the position)
	// a skip level of 3 allows to pass "CA" (for calls from identity )

	vector<string> tokens = MslTools::tokenize(_atomId, ",", true);
	if (tokens.size() == 1) {
		// no comma in string
		tokens = MslTools::tokenize( _atomId, " ");
		if (tokens.size() == 1) {
			// no comma in string
			tokens = MslTools::tokenize( _atomId, "_");
		}
	}
	_chainid = "";
	_resnum = 0;
	_icode = "";
	_identity = "";
	_atomName = "";
	if (tokens.size() != 4) {
		// not "A 37 ILE CA" format
		if (tokens.size() + _skiplevels < 4) {
			return false;
		} else {
			if (tokens.size() == 3 && _skiplevels >= 1) {
				// "37 ILE CA": add a blank chain
				tokens.insert(tokens.begin(), "");
			} else if (tokens.size() == 2 && _skiplevels >= 2) {
				// "ILE CA": add a blank chain a residue number
				tokens.insert(tokens.begin(), "0");
				tokens.insert(tokens.begin(), "");
			} else if (tokens.size() == 1 && _skiplevels >= 3) {
				// "CA": add a blank chain a residue number and identity
				tokens.insert(tokens.begin(), "");
				tokens.insert(tokens.begin(), "0");
				tokens.insert(tokens.begin(), "");
			} else {
				return false;
			}
		}
	}
	/*
	if (tokens.size() != 4) {
		_chainid = "";
		_resnum = 0;
		_icode = "";
		_identity = "";
		_atomName = "";
		_OK = false;
		return;
	}
	*/
	for (vector<string>::iterator k=tokens.begin(); k!=tokens.end(); k++) {
		*k = MslTools::trim(*k);
	}
	_chainid = tokens[0];
//	int resNum = 0;
//	string iCode;
	bool OK = MslTools::splitIntAndString(tokens[1], _resnum, _icode);
	_identity = tokens[2];
	_atomName = tokens[3];
	//_OK = true;
	return OK;
}

bool MslTools::compareAtomOfIdentityIds(string _id1, string _id2, unsigned int _skiplevels) {
	// check if two atoms ids are the same (even if they are separated by different
	// characters, it is not just a string compare, i.e. "A,37,ILE,CA" vs "A 37 ILE CA"
	string chainid1;
	int resnum1;
	string icode1;
	string identity1;
	string atomName1;
	string chainid2;
	int resnum2;
	string icode2;
	string identity2;
	string atomName2;
	if (!parseAtomOfIdentityId(_id1, chainid1, resnum1, icode1, identity1, atomName1, _skiplevels)) {
		return false;
	}
	if (!parseAtomOfIdentityId(_id2, chainid2, resnum2, icode2, identity2, atomName2, _skiplevels)) {
		return false;
	}
	if (chainid1 != chainid2 || resnum1 != resnum2 || icode1 != icode2 || identity1 != identity2 || atomName1 != atomName2) {
		return false;
	}
	return true;
}

// RegEx Functions  (only works with compile __BOOST__ flag on, otherwise returns false immediately)
bool MslTools::regex(string _lineToMatch, string _expression, vector<string> &_matches){

#ifndef __BOOST__
	return false;
#else
	boost::regex re(_expression);
	boost::cmatch matches;
	// dwkulp, note: doesn't regex_match, match the WHOLE string only?
	if (boost::regex_match(_lineToMatch.c_str(),matches,re)) {

	  // matches[0] contains the original string.  matches[n]
	  // contains a sub_match object for each matching
	  // subexpression
	  for (int i = 1; i < matches.size(); i++){

	    // sub_match::first and sub_match::second are iterators that
	    // refer to the first and one past the last chars of the
	    // matching subexpression
	    string match(matches[i].first, matches[i].second);	  

	    // Add match to list
	    _matches.push_back(match);
	  }

	} else{
	  return false;
	}

       // Successful matching...
        return true;
#endif
}
bool MslTools::regex(string _lineToMatch, string _expression, std::vector<std::pair<int,int> > &_matchIndicies){

#ifndef __BOOST__
	return false;
#else
	boost::regex re(_expression);

	boost::sregex_iterator r1(_lineToMatch.begin(),_lineToMatch.end(),re);
	boost::sregex_iterator r2;
	if (r1 == r2) return false;

	while (r1 != r2){
	  pair<int,int> a;
	  a.first = (*r1).position();
	  a.second = (*r1).position()+(*r1).length()-1;
	  _matchIndicies.push_back(a);

	  *r1++;
	}



       // Successful matching...
        return true;
#endif
}

std::string MslTools::stringf(const char * _format, ...){
  
  char buffer[1000];
  va_list args;
  va_start(args,_format);
  vsprintf(buffer,_format,args);
  va_end(args);


  return (string)buffer;
  
}
bool MslTools::replace(std::string &_string, const std::string &_replace, const std::string &_with, bool _replaceAll){

  size_t start_pos = _string.find(_replace);
  if(start_pos == std::string::npos)
    return false;


  start_pos = 0;
  while((start_pos = _string.find(_replace, start_pos)) != std::string::npos) {
    _string.replace(start_pos, _replace.length(), _with);
    start_pos += _with.length(); 
    if (!_replaceAll){
      break;
    }
  }

  return true;
}

double MslTools::getBoltzmannEnsembleEnergy(double _temp, std::vector<double>& _energies) {
	vector<double> probs = getBoltzmannProbabilities(_temp,_energies);
	double energy = 0;
	for(int i = 0; i < probs.size(); i++) {
		energy = energy + probs[i] * _energies[i];
	}
	return energy;
}

vector<double> MslTools::getBoltzmannProbabilities(double _temp, vector<double>& _energies) {
	vector<double> probs;
	double partition = 0;
	double minusRT = -MslTools::R * _temp;
	
	for(int i = 0; i < _energies.size(); i++) {
		probs.push_back(exp(_energies[i]/ minusRT));
		partition += probs.back();
	}

	for(int i = 0; i < probs.size(); i++) {
		probs[i] = probs[i]/partition;
	}
	return probs;
}

std::string MslTools::getBaseMSLDirectory() {
	const char* mslEnv = std::getenv("MSL_DIR");
	if (mslEnv == nullptr) {
		std::cerr << "Error: MSL_ENV environment variable is not set." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	return std::string(mslEnv);
}


