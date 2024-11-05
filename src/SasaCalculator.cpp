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
#include<SasaCalculator.h>

using namespace MSL;
using namespace std;


SasaCalculator::SasaCalculator() {
	setup(DEFAULT_OCCLUSION_POINTS, DEFAULT_PROBE_RADIUS);
}

SasaCalculator::SasaCalculator(AtomPointerVector& _atoms, double _probeRadius, const int _noOcclusionPoints) {
	deletePointers();
	setup(_noOcclusionPoints, _probeRadius);
	addAtoms(_atoms);	
}

vector<SasaAtom*> & SasaCalculator::getAtomPointers() {
	return sasaAtoms;
}

SasaCalculator::SasaCalculator(SasaCalculator &_sasaCalculator) {
	deletePointers();

	for(vector<SasaAtom*>::iterator i = _sasaCalculator.sasaAtoms.begin(); i != _sasaCalculator.sasaAtoms.end(); i++) {
		sasaAtoms.push_back( new SasaAtom(**i));
	}
	noOcclusionPoints = _sasaCalculator.noOcclusionPoints;
	atoms = _sasaCalculator.atoms;
	noOcclusionPoints = _sasaCalculator.noOcclusionPoints;
	useDefaultRadii = _sasaCalculator.useDefaultRadii;
	probeRadius = _sasaCalculator.probeRadius;
	atomRadii = _sasaCalculator.atomRadii;
	setTFactor = _sasaCalculator.setTFactor;
}


void SasaCalculator::setup(int _noOcclusionPoints, double _probeRadius) {
	useDefaultRadii = true;
	deletePointers();
	noOcclusionPoints = _noOcclusionPoints;
	probeRadius = _probeRadius;
	atomRadii.clear();
	atomRadii["C"] = 2.00;
	atomRadii["H"] = 1.09;
	atomRadii["N"] = 1.85;
	atomRadii["O"] = 1.735;
	atomRadii["CA"] = 1.71;
	atomRadii["FE"] = 0.65;
	atomRadii["S"] = 2.06;
	atomRadii["ZN"] = 1.09;
	atomRadii["HE"] = 1.48;
	atomRadii["NE"] = 1.53;
	atomRadii["NA"] = 1.36375;
	atomRadii["K"] = 1.76375;
	atomRadii["CL"] = 2.27;
	atomRadii["MG"] = 1.185;
	atomRadii["CS"] = 2.1;
	atomRadii["CU"] = 1.38;
	atomRadii["P"] = 2.15;
	setTFactor = false;
	
	// TODO : Find proper values for these
	//	SasaCalculator::atomRadii["F"] = 1.38;
	//	SasaCalculator::atomRadii["AS"] = 1.38;
		
//	AtomSasa.clear();
//	ResidueSasa.clear();
//	ChainSasa.clear();
}


void SasaCalculator::addAtoms(AtomPointerVector& _atoms) {
	atoms = _atoms;
	for(int i=0; i< atoms.size();i++) {
		if (useDefaultRadii) {
			if (atoms[i]->getRadius() <= 0.0) {
				if (atomRadii.find(atoms[i]->getElement()) != atomRadii.end()) {
					atoms[i]->setRadius(atomRadii[atoms[i]->getElement()]);
				}
			} else {
			}
		}
		sasaAtoms.push_back(new SasaAtom(atoms[i], probeRadius, 0));
	}
}

void SasaCalculator::deletePointers() {
	for(int i = 0; i < sasaAtoms.size(); i++) {
		delete sasaAtoms[i];	
	}
}


void SasaCalculator::calcSasa() {
	//vector<vector<SasaAtom*> > distList = findNeighbors(sasaAtoms);
	vector<AtomPointerVector> distList = findNeighbors();

//	AtomSasa.clear();
//	ResidueSasa.clear();
//	ChainSasa.clear();

	for(int curAtom = 0; curAtom < sasaAtoms.size(); curAtom++) {
		//Atom* thisAtom = sasaAtoms[curAtom]->getAtom();
		sasaAtoms[curAtom]->buildOcclusionPoints(noOcclusionPoints);

		bool purge_flag = false;
		for (unsigned int neighbour = 0; neighbour < distList[curAtom].size(); neighbour++) {
			if (neighbour == distList[curAtom].size() - 1) {
				purge_flag = true;
			}
			sasaAtoms[curAtom]->removeOccluded(*distList[curAtom][neighbour], purge_flag);	
			if (sasaAtoms[curAtom]->getSurfaceSphereCurrentSize() == 0) {
				break;
			}
		}

		atoms[curAtom]->setSasa(sasaAtoms[curAtom]->calcSasa());
		if (setTFactor) {
			atoms[curAtom]->setTempFactor(sasaAtoms[curAtom]->getSasa());
		}
	//	ChainSasa[atoms[curAtom]->getChainId()] += sasaAtoms[curAtom]->getSasa();
	//	string resnumIcode
	//	ResidueSasa[[atoms[curAtom]->getChainId()]
		sasaAtoms[curAtom]->deleteOcclusionPoints();		 // Delete the built occlusion points

	} 
}
 

string SasaCalculator::getSasaTable(bool _byAtom) {

	stringstream ss;
	if (_byAtom) {
		for (int i = 0 ; i < atoms.size();  i++) {
			char line[1000];
			//const Atom* A = sasaAtoms[i]->getAtom();
		//	sprintf(line,"%5d %3d %3s %4s %11.5f %11.5f %11.5f %1s %11.5f\n",i,atoms[i]->getResidueNumber(),atoms[i]->getResidueName().c_str(),atoms[i]->getName().c_str(),atoms[i]->getX(),atoms[i]->getY(),atoms[i]->getZ(),atoms[i]->getChainId().c_str(),sasaAtoms[i]->getSasa());
			sprintf(line,"%5d %1s %4d%1s %3s %-4s %11.5f\n",i+1,atoms[i]->getChainId().c_str(),atoms[i]->getResidueNumber(),atoms[i]->getResidueIcode().c_str(),atoms[i]->getResidueName().c_str(),atoms[i]->getName().c_str(),sasaAtoms[i]->getSasa());
			ss << line;
		}
	} else {
	       /*
	SASA reference:
	Protein Engineering vol.15 no.8 pp.659–667, 2002
	Quantifying the accessible surface area of protein residues in their local environment
	Uttamkumar Samanta Ranjit P.Bahadur and  Pinak Chakrabarti
      */
      map<string,double> refSasa;
      refSasa["G"] = 83.91;
      refSasa["A"] = 116.40;
      refSasa["S"] = 125.68;
      refSasa["C"] = 141.48;
      refSasa["P"] = 144.80;
      refSasa["T"] = 148.06;
      refSasa["D"] = 155.37;
      refSasa["V"] = 162.24;
      refSasa["N"] = 168.87;
      refSasa["E"] = 187.16;
      refSasa["Q"] = 189.17;
      refSasa["I"] = 189.95;
      refSasa["L"] = 197.99;
      refSasa["H"] = 198.51;
      refSasa["K"] = 207.49;
      refSasa["M"] = 210.55;
      refSasa["F"] = 223.29;
      refSasa["Y"] = 238.30;
      refSasa["R"] = 249.26;
      refSasa["W"] = 265.42;

                
		map<string,double> residueSasa;
		string resId = "";
		char line[1000];

		for (int i = 0 ; i < atoms.size();  i++) {
			sprintf(line,"      %1s %4d%1s %3s      ",atoms[i]->getChainId().c_str(),atoms[i]->getResidueNumber(),atoms[i]->getResidueIcode().c_str(),atoms[i]->getResidueName().c_str());
			resId = string(line,23);
			//cout << resId << endl;
			if(residueSasa.find(resId) != residueSasa.end()) {
				residueSasa[resId] += sasaAtoms[i]->getSasa();
			} else {
				residueSasa[resId] = sasaAtoms[i]->getSasa();
			}
		}

		for(map<string,double>::iterator it = residueSasa.begin(); it != residueSasa.end(); it++) {
		  string chain;
		  int resnum;
		  string icode;
		  string identity;
		  MslTools::parseIdentityId(it->first, chain,resnum,icode,identity);
		  map<string,double>::iterator findIt = refSasa.find(MslTools::getOneLetterCode(identity));
		  double nSasa = 0;
		  if (findIt != refSasa.end()){
		    nSasa = it->second / refSasa[MslTools::getOneLetterCode(identity)];
		  } else {
		    nSasa = 1;
		  }
			sprintf (line,"%11.5f\t%11.5f\n",it->second,nSasa);
			ss << it->first << line;
		}

		/*
		string prevChain = "";
		int prevRes = 0;
		string prevIcode = "";
		double sasa = 0.0;
		for (int i = 0 ; i < atoms.size();  i++) {
			string chain = atoms[i]->getChainId();
			int res = atoms[i]->getResidueNumber();
			string icode = atoms[i]->getResidueIcode();
			if (res != prevRes || chain != prevChain || icode != prevIcode) {
				if (i != 0) {
					char line[1000];
					sprintf(line,"      %1s %4d%1s %3s      %11.5f\n",atoms[i-1]->getChainId().c_str(),atoms[i-1]->getResidueNumber(),atoms[i-1]->getResidueIcode().c_str(),atoms[i-1]->getResidueName().c_str(),sasa);
					ss << line;
				}
				sasa = 0.0;
			}
			sasa += sasaAtoms[i]->getSasa();
			if (i == atoms.size() - 1) {
				char line[1000];
				sprintf(line,"      %1s %4d%1s %3s      %11.5f\n",atoms[i]->getChainId().c_str(),atoms[i]->getResidueNumber(),atoms[i]->getResidueIcode().c_str(),atoms[i]->getResidueName().c_str(),sasa);
				ss << line;
			}
			prevChain = chain;
			prevRes = res;
			prevIcode = icode;
		}
		*/
	}
	return ss.str();
}

vector<AtomPointerVector> SasaCalculator::findNeighbors() {
	vector<AtomPointerVector> distList(atoms.size());
	//vector<vector<SasaAtom*> > distList(_sasaAtoms.size());
	vector<vector<double> > distVals(atoms.size());
//	unsigned int dCounter = 0;
	if (atoms.size() < DEFAULT_MIN_ATOMSIZE_FOR_ATOMGRID) {
		for (int i = 0; i< atoms.size(); i++) {
			double r1 = atoms[i]->getRadius();
			if (r1 <= 0) {
				continue;
			}

			// build a list of neighbors for each atom
			for (unsigned int j=i+1; j < atoms.size(); j++) {

				double r2 = atoms[j]->getRadius();
				if (r2 <= 0.0) {
					continue;
				}
				double d = atoms[i]->distance(*atoms[j]);
		//		dCounter++;
				if(d < (r1 + r2 + 2 * probeRadius )) {
					distList[i].push_back(atoms[j]);
					distVals[i].push_back(d);
					distList[j].push_back(atoms[i]);
					distVals[j].push_back(d);
				} 
			}
			// sort the lists by distance, this way the majority of the surface is
			// taken away immediately, speeding up the process
			vector<unsigned int> index(distVals[i].size(), 0);
			for (unsigned int j=0; j<index.size(); j++) {
				index[j] = j;
			}
			MslTools::quickSortWithIndex(distVals[i], index);
			AtomPointerVector tmpDistList (distList[i].size(), NULL);
			for (unsigned int j=0; j<index.size(); j++) {
				tmpDistList[j] = distList[i][ index[j] ];
			}
			distList[i] = tmpDistList;
		}
	} else {
		/*
		double rMax1;
		double rMax2;
		for (int i = 0; i< atoms.size(); i++) {
			double r = atoms[i]->getRadius();
			if (i==0 || r >= rMax1) {
				rMax2 = rMax1;
				rMax1 = r;
			} else if (i==1 || r >= rMax2) {
				rMax2 = r;
			}
		}
		cout << "UUUU " << rMax1 << " " << rMax2 << endl;
		Atom3DGrid grid(atoms, rMax1+rMax2+2*probeRadius);
		*/
		Atom3DGrid grid(atoms, 5+2*probeRadius);
		for (int i = 0; i< atoms.size(); i++) {
			double r1 = atoms[i]->getRadius();
			if (r1 <= 0) {
				continue;
			}
			AtomPointerVector neighbors = grid.getNeighbors(i);	
			// build a list of neighbors for each atom
			for (unsigned int j=0; j < neighbors.size(); j++) {

				double r2 = neighbors[j]->getRadius();
				if (r2 <= 0.0) {
					continue;
				}
				double d = atoms[i]->distance(*neighbors[j]);
			//	dCounter++;
				if(d < (r1 + r2 + 2 * probeRadius )) {
					distList[i].push_back(neighbors[j]);
					distVals[i].push_back(d);
					//distList[j].push_back(atoms[i]);
					//distVals[j].push_back(d);
				} 
			}
			// sort the lists by distance, this way the majority of the surface is
			// taken away immediately, speeding up the process
			vector<unsigned int> index(distVals[i].size(), 0);
			for (unsigned int j=0; j<index.size(); j++) {
				index[j] = j;
			}
			MslTools::quickSortWithIndex(distVals[i], index);
			AtomPointerVector tmpDistList (distList[i].size(), NULL);
			for (unsigned int j=0; j<index.size(); j++) {
				tmpDistList[j] = distList[i][ index[j] ];
			}
			distList[i] = tmpDistList;

		}
	}
	//cout << "UUUU d counter " << dCounter << endl;
	return distList;
}

void SasaCalculator::printSasaTable(bool _byAtom) {
	cout << getSasaTable(_byAtom);
}


double SasaCalculator::getTotalSasa(){
	double sasa = 0.0;
	for (int i = 0 ; i < atoms.size();  i++) {
		sasa += sasaAtoms[i]->getSasa();
	}

	return sasa;
}
double SasaCalculator::getResidueSasa(string _positionId){

  double sasa = 0.0;
  for (int i = 0 ; i < atoms.size();  i++) {
    string atomPositionId = MslTools::getPositionId(atoms[i]->getChainId(),atoms[i]->getResidueNumber(),atoms[i]->getResidueIcode());
    if (_positionId == atomPositionId){
      sasa += sasaAtoms[i]->getSasa();
    }
  }

  return sasa;
}
