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

#include "CIFFormat.h"

using namespace MSL;
using namespace std;


/*

Here is the mmCIF format for the _atom_site category (there are more fields on-line):

loop_
_atom_site.group_PDB 
_atom_site.id 
_atom_site.type_symbol 
_atom_site.label_atom_id 
_atom_site.label_alt_id 
_atom_site.label_comp_id 
_atom_site.label_asym_id 
_atom_site.label_entity_id 
_atom_site.label_seq_id 
_atom_site.pdbx_PDB_ins_code 
_atom_site.Cartn_x 
_atom_site.Cartn_y 
_atom_site.Cartn_z 
_atom_site.occupancy 
_atom_site.B_iso_or_equiv 
_atom_site.pdbx_formal_charge 
_atom_site.auth_seq_id 
_atom_site.auth_comp_id 
_atom_site.auth_asym_id 
_atom_site.auth_atom_id 
_atom_site.pdbx_PDB_model_num 
ATOM   1    N N   . TYR A 1 9   ? -16.253 28.841 -23.028 1.00 36.69  ? 9   TYR A N   1 
ATOM   2    C CA  . TYR A 1 9   ? -16.407 29.570 -21.787 1.00 44.83  ? 9   TYR A CA  1 
ATOM   3    C C   . TYR A 1 9   ? -16.522 28.610 -20.612 1.00 43.19  ? 9   TYR A C   1 
ATOM   4    O O   . TYR A 1 9   ? -16.008 28.900 -19.551 1.00 46.80  ? 9   TYR A O   1 
 */

Atom * CIFFormat::createAtomFromAtomSiteLine(const string &_pdbAtomLine, vector<string> &_fields, bool _allowPartialRead, bool _useAuthFields){
  
  vector<string> toks = MslTools::tokenizeAndTrim(_pdbAtomLine," ");

  Atom *at = new Atom();

  bool _fillAuthAtomName = false;
  bool _fillAuthResName  = false;
  bool _fillAuthChainId  = false;
  bool _fillAuthResNum   = false;
  for (uint i = 0; i < toks.size();i++){

    //cout << "FIELD: "<<_fields[i]<<" toks: "<<toks[i]<<endl;
    
    // Skip undef values, try to use "auth" values, if not then it will keep default values from Atom constructor

    if (toks[i] == "." || toks[i] == "?"){

      if (_fields[i] == "label_atom_id") _fillAuthAtomName = true;
      if (_fields[i] == "label_comp_id") _fillAuthResName  = true;
      if (_fields[i] == "label_asym_id") _fillAuthChainId  = true;
      if (_fields[i] == "label_seq_id")  _fillAuthResNum   = true;
      continue;
      
    } else {
    
      try {
	
	if (_fields[i] == "type_symbol")                       at->setElement(toks[i]);
	if (_fields[i] == "label_atom_id")                     at->setName(toks[i]);
	if (_fields[i] == "label_comp_id")                     at->setResidueName(toks[i]);
	if (_fields[i] == "label_asym_id")                     at->setChainId(toks[i]);
	if (_fields[i] == "label_seq_id")                      at->setResidueNumber(MslTools::toInt(toks[i], "Residue Number is not an int"));
	if (_fields[i] == "pdbx_PDB_ins_code")                 at->setResidueIcode(toks[i]);
	if (_fields[i] == "Cartn_x")                           at->getAllCoor()[0]->setX(MslTools::toDouble(toks[i], "X is not a double"));
	if (_fields[i] == "Cartn_y")                           at->getAllCoor()[0]->setY(MslTools::toDouble(toks[i], "Y is not a double"));
	if (_fields[i] == "Cartn_z")                           at->getAllCoor()[0]->setZ(MslTools::toDouble(toks[i], "Z is not a double"));
	if (_fields[i] == "B_iso_or_equiv")                    at->setTempFactor(MslTools::toDouble(toks[i], "TempFactor is not a double"));


	if (_fields[i] == "auth_atom_id" && (_fillAuthAtomName || _useAuthFields)) at->setName(toks[i]);
	if (_fields[i] == "auth_comp_id" && (_fillAuthResName  || _useAuthFields)) at->setResidueName(toks[i]);
	if (_fields[i] == "auth_asym_id" && (_fillAuthChainId  || _useAuthFields)) at->setChainId(toks[i]);
	if (_fields[i] == "auth_seq_id"  && (_fillAuthResNum   || _useAuthFields)) at->setResidueNumber(MslTools::toInt(toks[i], "Residue Number2 is not an int"));
	
	// SegId is in different field.
	
      } catch(exception &e){
	cerr << "ERROR 34919 CIFFormat createAtomFromAtomSiteLine "<<e.what()<<" line is: "<<_pdbAtomLine<<endl;
	exit(34919);
      }
      

      
    }

  }

  return at;

}


string CIFFormat::createLoopAtomSite(AtomPointerVector &_ats, vector<string> &_fields, bool _addHeader, bool _addTail){

  std::stringstream a;

  if (_fields.size() == 0){
    _fields.push_back("_atom_site.group_PDB");
    _fields.push_back("_atom_site.id");
    _fields.push_back("_atom_site.type_symbol"); 
    _fields.push_back("_atom_site.label_atom_id"); 
    _fields.push_back("_atom_site.label_comp_id");  
    _fields.push_back("_atom_site.label_asym_id");   
    _fields.push_back("_atom_site.label_seq_id");    
    _fields.push_back("_atom_site.pdbx_PDB_ins_code"); 
    _fields.push_back("_atom_site.Cartn_x"); 
    _fields.push_back("_atom_site.Cartn_y"); 
    _fields.push_back("_atom_site.Cartn_z");
    _fields.push_back("_atom_site.occupancy");
    _fields.push_back("_atom_site.B_iso_or_equiv"); 
  }
  if (_addHeader){
    a << "loop_\n";
    for (uint i = 0; i < _fields.size();i++){
      a << _fields[i] << "\n";
    }
  }
  
  for (uint i = 0; i < _ats.size();i++){
    a << createAtomSiteLineFromAtom(_ats(i), _fields,i);
  }

  if (_addTail){
    a << "#\n";
  }

  return a.str();
}

string CIFFormat::createAtomSiteLineFromAtom(Atom &_at, vector<string> &_fields, unsigned int _atomNumber){

  unsigned int atomNumber = _atomNumber;
  if (_atomNumber == 0) atomNumber = 1;
  
  if (_fields.size() == 0){
    _fields.push_back("_atom_site.group_PDB");
    _fields.push_back("_atom_site.id");
    _fields.push_back("_atom_site.type_symbol"); 
    _fields.push_back("_atom_site.label_atom_id"); 
    _fields.push_back("_atom_site.label_comp_id");  
    _fields.push_back("_atom_site.label_asym_id");   
    _fields.push_back("_atom_site.label_seq_id");    
    _fields.push_back("_atom_site.pdbx_PDB_ins_code"); 
    _fields.push_back("_atom_site.Cartn_x"); 
    _fields.push_back("_atom_site.Cartn_y"); 
    _fields.push_back("_atom_site.Cartn_z");
    _fields.push_back("_atom_site.occupancy");
    _fields.push_back("_atom_site.B_iso_or_equiv"); 
  }


  std::stringstream a;
  a << "ATOM ";
  a << atomNumber << " ";
  for (uint i = 0; i < _fields.size();i++){
    if (_fields[i] == "_atom_site.type_symbol")             { if (_at.getElement() == "") { a << "? "; } else { a << _at.getElement() << " ";} }
    if (_fields[i] == "_atom_site.label_atom_id")	    { if (_at.getName() == "") { a << "? "; } else {  a << _at.getName() << " ";} }
    if (_fields[i] == "_atom_site.label_comp_id")	    { if (_at.getResidueName() == "") { a << "? "; } else { a << _at.getResidueName() << " ";} }
    if (_fields[i] == "_atom_site.label_asym_id")	    { if (_at.getChainId() == "") { a << "? "; } else { a << _at.getChainId() << " ";} }
    if (_fields[i] == "_atom_site.label_seq_id")	    { a << _at.getResidueNumber() << " ";}
    if (_fields[i] == "_atom_site.pdbx_PDB_ins_code")       { if (_at.getResidueIcode() == "") { a << "? "; } else { a << _at.getResidueIcode() << " "; } }
    if (_fields[i] == "_atom_site.Cartn_x")	            { a << _at.getX() << " "; } 
    if (_fields[i] == "_atom_site.Cartn_y")	            { a << _at.getY() << " ";} 
    if (_fields[i] == "_atom_site.Cartn_z")	            { a << _at.getZ() << " ";}
    if (_fields[i] == "_atom_site.occupancy")	            { a << 1.00 << " "; }
    if (_fields[i] == "_atom_site.B_iso_or_equiv")	    { a << _at.getTempFactor() << " ";}
  }

  a << "\n";

  
  return a.str();
  
}

std::string CIFFormat::createLoopAtomSite(AtomPointerVector &_ats, bool _addHeader, bool _addTail){
  std::vector<std::string> _fields;
  return createLoopAtomSite(_ats, _fields, _addHeader, _addTail);
}

std::string CIFFormat::createAtomSiteLineFromAtom(Atom &_at,unsigned int _atomNumber){
  std::vector<std::string> _fields;
  return createAtomSiteLineFromAtom(_at, _fields, _atomNumber);
}

