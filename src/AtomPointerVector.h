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


#ifndef ATOMPOINTERVECTOR_H
#define ATOMPOINTERVECTOR_H

// STL Includes
#include <vector>
#include <sstream>
#include <sys/types.h>

// MSL Includes
#include "Atom.h"
#include "CartesianPoint.h"

// BOOST Includes
#ifdef __BOOST__
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#endif

// Forward Declarations
namespace MSL { 
//class Atom;


class AtomPointerVector : public std::vector<Atom *> {
	public:
		AtomPointerVector();
		AtomPointerVector(unsigned int _size, Atom * _pointer=NULL);
		AtomPointerVector(const AtomPointerVector & _atoms);
		~AtomPointerVector();

		void operator=(const AtomPointerVector & _atoms);
		AtomPointerVector operator+(const AtomPointerVector & _atoms);
		AtomPointerVector operator-(const AtomPointerVector & _atoms);
		void operator+=(const AtomPointerVector & _atoms);
		void operator-=(const AtomPointerVector & _atoms);

		Atom & operator()(unsigned int _n);

		friend std::ostream & operator<<(std::ostream &_os, const AtomPointerVector & _av)  {_os << _av.toString(); return _os;};

		// Setter/Getter Functions
		void   setName(std::string _name);
		std::string getName() const;
	//	CartesianPoint getGeometricCenter() const;
		CartesianPoint& getGeometricCenter(unsigned int _stamp=0);

		// Geometric Center
	 //       void updateGeometricCenter(unsigned int _updateStamp=0);

		double rmsd(const AtomPointerVector &_av) const;

        int clashCheck(AtomPointerVector &_av, double _dist2_tol=9);

		// Get the minimum or maximum number of alt coordinates for over all atoms
		int getMinAltConf();
		int getMaxAltConf();

		void setActiveConf(int _conf);

        void setTempFactor(double &_temp);
	     //   void translate(double _x, double _y, double _z);         
	     //   void translate(const CartesianPoint &_vec);         
		//void rotate(const Matrix &_rotMat);


		/***************************************************
		 *  Saving coordinates to buffers:
		 *
		 *  coordinates can be saved to named buffers (string _coordName),
		 *  and copied back from them
		 *
		 *  The difference between save coordinates to a buffer, and 
		 *  having multiple alternate coor is that the saved coord 
		 *  are simply a buffer that can be restored
		 *
		 *  Coor can be saved to buffer with two different commands:
		 *    saveCoor:
		 *      - saveCoor saves ONLY the current coor
		 *      - when restored with applySavedCoor, a buffer created with
		 *        saveCoor will replace the CURRENT coorinate only
		 *    saveAltCoor:
		 *      - saveAltCoor saves ALL alternative coordinates and
		 *        also remembers what was the current coordinate
		 *      - when restored with the same applySavedCoor, a buffer
		 *        created with saveAltCoor will wipe off all alternative
		 *        cordinates and recreate the situation that was present
		 *        when the buffer was saved
		 *
		 *  More details in Atom.h
		 ***************************************************/
		void saveCoor(std::string _coordName);
		void saveAltCoor(std::string _coordName);
		bool applySavedCoor(std::string _coordName);
		void clearSavedCoor(std::string _coordName="");		

		//void addAtomRanking(std::string _atomKey, double _val);
		//void sort();

		std::string toString() const;

		void deletePointers();


		// Boost related, is ok if no BOOST libraries are being used, just a std::string.
		void setArchiveType(std::string _type) { archiveType = _type; }
		std::string getArchiveType() { return archiveType; }


		std::vector< std::vector< std::map<std::string, Atom*> > > subdivideByChainAndPosition();
		std::vector< std::vector< std::map<std::string, std::map<std::string, Atom*> > > > subdivideByChainPositionAndIndentity();

	private:
		void setup();

		std::string name;
		CartesianPoint geometricCenter;	

//		std::map<std::string, double> atomRankings;
		unsigned int updateStamp;

		std::string archiveType;



		// BOOST-RELATED FUNCTIONS , keep them away from main class def.
#ifdef __BOOST__
	public:

		void save_checkpoint(std::string filename) const{

			if (archiveType == "binary"){
				std::ofstream fout(filename.c_str(),std::ios::binary);
				boost::archive::binary_oarchive oa(fout);
				oa << (*this);
			} else if (archiveType == "xml"){
				std::ofstream fout(filename.c_str());
				boost::archive::xml_oarchive oa(fout);
				oa << boost::serialization::make_nvp("AtomPointerVector",*this);
			} else {
				std::ofstream fout(filename.c_str());
				boost::archive::text_oarchive oa(fout);
				oa << (*this);
			}

		}

		void load_checkpoint(std::string filename){

			if (archiveType == "binary"){
				std::ifstream fin(filename.c_str(), std::ios::binary);
				boost::archive::binary_iarchive ia(fin);
				ia >> (*this);
			} else if (archiveType == "xml"){
				std::ifstream fin(filename.c_str());
				boost::archive::xml_iarchive ia(fin);
				ia >> boost::serialization::make_nvp("AtomPointerVector",*this);
			} else {
				std::ifstream fin(filename.c_str());
				boost::archive::text_iarchive ia(fin);
				ia >> (*this);
			}
		}
	private:

		friend class boost::serialization::access;		


		template<class Archive> void serialize(Archive & ar, const unsigned int version){

			using boost::serialization::make_nvp;
			//ar & boost::serialization::base_object<std::vector<Atom *> >(*this);
			//ar & geometricCenter;

			ar & make_nvp("name",name);
			ar & make_nvp("atoms",boost::serialization::base_object<std::vector<Atom *> >(*this));
		}
#else
	public:
		void save_checkpoint(std::string filename) const{
			std::cout << "NO IMPLEMENTATION OF SAVE_CHECKPOINT WITHOUT BOOST LIBRARIES INSTALLED.\n";
		}
		void load_checkpoint(std::string filename) const{
			std::cout << "NO IMPLEMENTATION OF LOAD_CHECKPOINT WITHOUT BOOST LIBRARIES INSTALLED.\n";
		}
#endif
};

// INLINE FUNCTIONS
inline Atom & AtomPointerVector::operator()(unsigned int _n) { return *((*this)[_n]); }
inline void AtomPointerVector::setName(std::string _name) { name = _name; }
inline std::string AtomPointerVector::getName() const { return name;  }
//inline CartesianPoint AtomPointerVector::getGeometricCenter() const { return geometricCenter; }
inline std::string AtomPointerVector::toString() const {
	std::stringstream ss;
	for (AtomPointerVector::const_iterator k = begin(); k!=end(); k++) {
		ss << **k << std::endl;
	}
	/*
	for (uint i = 0; i < (*this).size();i++){
		ss << "I: "<<i<<std::endl;
		ss << (*this)(i)<<std::endl;
	}
	*/
	return ss.str();
}

inline void AtomPointerVector::setTempFactor(double &_temp){
  //for (AtomPointerVector::iterator k = begin(); k!=end(); k++) { (*k).setTempFactor(_temp);}
  for (uint i = 0; i < size();i++){ (*this)(i).setTempFactor(_temp);}
}
  
}

#endif
