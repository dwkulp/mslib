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

#include "PSSMCreator.h"
#include <cctype>

using namespace MSL;
using namespace std;

#include "MslOut.h"
static MslOut MSLOUT("PSSMCreator");

PSSMCreator::PSSMCreator(){

  expectedValues["A"] = 1.0;
  expectedValues["C"] = 1.0;
  expectedValues["D"] = 1.0;
  expectedValues["E"] = 1.0;
  expectedValues["F"] = 1.0;
  expectedValues["G"] = 1.0;
  expectedValues["H"] = 1.0;
  expectedValues["I"] = 1.0;
  expectedValues["K"] = 1.0;
  expectedValues["L"] = 1.0;
  expectedValues["M"] = 1.0;
  expectedValues["N"] = 1.0;
  expectedValues["P"] = 1.0;
  expectedValues["Q"] = 1.0;
  expectedValues["R"] = 1.0;
  expectedValues["S"] = 1.0;
  expectedValues["T"] = 1.0;
  expectedValues["U"] = 1.0;
  expectedValues["V"] = 1.0;
  expectedValues["W"] = 1.0;
  expectedValues["Y"] = 1.0;
  expectedValues["X"] = 1.0;
}

PSSMCreator::~PSSMCreator(){}


void PSSMCreator::addMultipleSequenceAlignment(string _fastaFile, string _regex){
  fastaFile = _fastaFile;
  fin.open(fastaFile);
  fin.read(_regex);
  fin.close();
  sequences = fin.getSequences();
}
void PSSMCreator::setSequences(map<string,string> &_seqs){
  sequences = _seqs;
}

void PSSMCreator::create(PSSMType _type){
  if (sequences.size() == 0){
    cerr << "ERROR 9732 no sequences in PSSMCreator"<<endl;
    exit(9732);
  }
  map<string,string>::iterator it;

  int len = sequences.begin()->second.length();
  observedValues.clear();
  observedValues.resize(len);
  vector<double> totalAA;
  totalAA.resize(len);
  for (it = sequences.begin();it != sequences.end();it++){
    if (it->second.length() != len){
      cerr << "ERROR sequence with name "<<it->first<< " has length "<<it->second.length()<< " and the first sequence with name "<<sequences.begin()->first<<" has length "<<sequences.begin()->second.length()<<" skipping this sequence and not adding its values to the PSSM"<<endl;
      continue;
    }
    for (uint i = 0; i < it->second.length();i++){
      string aa = MslTools::stringf("%c",it->second[i]);
      observedValues[i][aa] += 1.0;
      if (aa != "-"){
	totalAA[i] += 1.0;
      }
    }
  }

  // Now create log-odds PSSM
  scoreFunction.clear();
  scoreFunction.resize(observedValues.size());
  map<string,double>::iterator it1;
  for (uint i = 0; i < observedValues.size();i++){
    scoreFunction[i]["entropy"] = 0.0;
    scoreFunction[i]["mostfreq"] = 0.0;

    for (it1 = observedValues[i].begin();it1 != observedValues[i].end();it1++){

      map<string,double>::iterator it2;
      it2 = expectedValues.find(it1->first);
      if (it2 == expectedValues.end()){
	cerr << "No expected value for "<<it1->first<<endl;
      } else {
	switch (_type) {
	   case logodds:
	     scoreFunction[i][it1->first] = -log( (it1->second/totalAA[i])/expectedValues[it1->first]);	     
	     break;
	   case freq:
	     scoreFunction[i][it1->first] = it1->second/totalAA[i];
	     break;
	}
	
	// Compute Shannon Entropy for this position 'i'
	scoreFunction[i]["entropy"] += ( (it1->second/totalAA[i]) * log( (it1->second/totalAA[i]) ) ); 
	MSLOUT.stream() << "Position ["<<i<<"]: "<< it1->first<<" "<<it1->second<<" "<<totalAA[i]<<" "<<(it1->second/totalAA[i])<<" "<<(it1->second/totalAA[i])*log( (it1->second/totalAA[i]))<<" "<<scoreFunction[i]["entropy"]<<endl;
					 

	// Keep track of the most frequent AA at position 'i'
	if (scoreFunction[i]["mostfreq"] < it1->second/totalAA[i]){
	  scoreFunction[i]["mostfreq"]  = it1->second/totalAA[i];
	}

      }
    }

    scoreFunction[i]["entropy"] *= -1;
    MSLOUT.stream() << "ENTROPY["<<i<<"] = "<<scoreFunction[i]["entropy"]<<endl;
  }


}

vector<double>  PSSMCreator::getScoreFunction(string _sequence, string _nameRefSeq, int _beginOffset, string _specificScoreType){


  string newSequence = "";
  vector<double> results;
  if (_nameRefSeq != ""){

    // Find Reference sequence
    string refSeq = getSequence(_nameRefSeq);
    if (refSeq == ""){
      cerr << "ERROR PSSMCreator::getScoreFunction() couldn't find name: "<<_nameRefSeq<<" in fasta file"<<endl;
      return results;
    }

    // Create aligned sequence (add non-alpha characters to input _sequence)
    int index = 0;
    for (uint i = 0; i < refSeq.size();i++){

      if (i < _beginOffset) {
	newSequence += "-";
	continue;
      }

      if (isalpha(refSeq[i])){
	if (index >= _sequence.length()) {
	  newSequence += "-";
	} else {
	  newSequence += MslTools::stringf("%c",_sequence[index++]);
	}
      } else {
	newSequence += MslTools::stringf("%c",refSeq[i]);
      }

    } // END refSeq

    MSLOUT.stream() << "Ref     Sequence: "<<refSeq<<endl;
    MSLOUT.stream() << "Aligned Sequence: "<<newSequence<<endl;
  } // END refSeq != ""





  if (newSequence.size() != scoreFunction.size()){
    cerr << "ERROR input sequence and score function have different lengths: "<<newSequence.size()<<" "<<scoreFunction.size()<<endl;
    return results;
  }
  for (uint i = 0; i < newSequence.size();i++){

    // Skip dashes and * and whatnot...
    if (isalpha(newSequence[i])){
      string scoreType = MslTools::stringf("%c",newSequence[i]);
      if (_specificScoreType != ""){
	scoreType = _specificScoreType;
      }

      results.push_back(scoreFunction[i][scoreType]);
      fprintf(stdout,"score of %d is %8.3f amino acid is: %c at index %d\n",i,scoreFunction[i][scoreType],newSequence[i],int(results.size()));
    }

  }
  
  if (results.size() != _sequence.size()){
    cerr << "ERROR 1834 PSSMCreator::getScoreFunction() result vector has length: "<<results.size()<<" but input sequence has length: "<<_sequence.size()<<endl;
  }

  

  return results;
}

void PSSMCreator::readReferenceCounts(string _filename){
  
  referenceFile = _filename;
  map<string,double> refCounts;
  vector<string> lines;
  if (!MslTools::readTextFile(lines,referenceFile)){
    return;
  }

  for (uint i = 0; i < lines.size();i++){
    if (lines[i].size() == 0) continue;// skip blank lines
    vector<string> toks = MslTools::tokenize(lines[i]," ");
    if (toks.size() != 2){
      cerr << "ERROR 341234 PSSMCreator::readReferenceCounts() line of refCount file: "<<_filename<<" has more than 2 whitespace-delmited tokens: "<<toks.size()<<" line is: "<<lines[i]<<endl;
      continue;
    }
    refCounts[toks[0]] = MslTools::toDouble(toks[1]);
    
  }

  addReferenceCounts(refCounts);
}
void PSSMCreator::addReferenceCounts(map<string,double> &_referenceCounts){
  
  map<string,double>::iterator it;
  double sum = 0.0;
  for (it = _referenceCounts.begin();it != _referenceCounts.end();it++){
    sum += it->second;
  }

  for (it = _referenceCounts.begin();it != _referenceCounts.end();it++){
    expectedValues[it->first] = (it->second / sum);
  }

}
vector<map<string,double> > PSSMCreator::getFrequencies(){

  vector<map<string,double> > results;
  for (uint i = 0; i < scoreFunction.size();i++){
      results.push_back(scoreFunction[i]);      
  }

  return results;
  
}
vector<map<string,double> > PSSMCreator::getFrequencies(string _nameRefSeq, int _resiBegin, int _resiEnd){

  if (_resiBegin < 0 || _resiBegin >= scoreFunction.size()) {
    cerr << "ERROR 9234 PSSMCreator::getFrequenices() _resiBegin is out of range: "<<_resiBegin<<" size of scoreFunction: "<<scoreFunction.size()<<endl;
    exit(9234);
  }
  if (_resiEnd < 0 || _resiEnd >= scoreFunction.size()) {
    cerr << "ERROR 9234 PSSMCreator::getFrequenices() _resiEnd is out of range: "<<_resiEnd<<" size of scoreFunction: "<<scoreFunction.size()<<endl;
    exit(9234);
  }

  // Figure out correct begin/end indicies in internal numbering
  int actualResiBegin = getInternalNumbering(_nameRefSeq, _resiBegin);
  int actualResiEnd   = getInternalNumbering(_nameRefSeq, _resiEnd);

  // Find Reference sequence
  string refSeq = getSequence(_nameRefSeq);

  MSLOUT.stream() << "Actual range: "<<actualResiBegin<<"-"<<actualResiEnd<< " "<<scoreFunction.size()<<endl;
  vector<map<string,double> > results;
  for (uint i = actualResiBegin; i <= actualResiEnd;i++){

    if (isalpha(refSeq[i])){
      results.push_back(scoreFunction[i]);      
      MSLOUT.stream() << "ScoreFunction at "<<i<<" has first AA["<<scoreFunction[i].begin()->first<< "]"<<endl;
      MSLOUT.stream() << "Freq: "<<scoreFunction[i].begin()->second<<endl;
    }

  }

  return results;
  
}

string PSSMCreator::getSequence(string _nameRefSeq, int _resiStart, int _resiEnd){

  int actualStart = getInternalNumbering(_nameRefSeq,_resiStart);
  int actualEnd   = getInternalNumbering(_nameRefSeq,_resiEnd);

  // Find Reference sequence
  string refSeq = getSequence(_nameRefSeq);
  if (actualEnd > refSeq.size()){
    cerr << "ERROR 9933 PSSMCreator::getSequence() actualEnd is greater than ref seq size: "<<refSeq.size()<<" actual end: "<<actualEnd<<" input end: "<<_resiEnd<<endl;
    exit(9933);
  }

  string result = "";
  for (uint i = actualStart; i <= actualEnd;i++){
    if (isalpha(refSeq[i])){
      result += refSeq[i];
    }
  }

  return result;
}
int PSSMCreator::getInternalNumbering(string _nameRefSeq, int _resi){
  
  // Find Reference sequence
  string refSeq = getSequence(_nameRefSeq);
  if (refSeq == ""){
    cerr << "ERROR 1231 PSSMCreator::getInternalNumbering() couldn't find name: "<<_nameRefSeq<<" in fasta file"<<endl;
    exit(1231);
  }

  MSLOUT.stream() << "Get Internal Numbering: "<<_resi<<" from "<<_nameRefSeq<<" size "<<refSeq.size()<<endl;
  int actualResi = 0;
  int i = 0;
  do {

    // Increment counter when character is present 
    if (isalpha(refSeq[actualResi++])){
      i++;
    }

  } while (i <  _resi);
  
  actualResi--;

  return actualResi;
}


string PSSMCreator::getSequence(string _name){
  map<string,string>::iterator it;
  it = sequences.find(_name);
  if (it == sequences.end()){
    return "";
  }

  return it->second;
}
