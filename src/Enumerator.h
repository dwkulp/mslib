/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Simulation Library)n
 Copyright (C) 2009 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan

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

#ifndef ENUMERATOR_H
#define ENUMERATOR_H

#include <vector>
#include <iostream>
#include <stdio.h>
#include <cstdlib>



/*! \brief This class is used to enumerate combinations of states
 */


namespace MSL { 
class Enumerator {
	public:
		Enumerator();
		Enumerator(std::vector<unsigned int> _states);
		Enumerator(std::vector<std::vector<unsigned int> > & _values);
		Enumerator(Enumerator & _enum);
		~Enumerator();

		std::vector<unsigned int> & operator[](size_t n);
		void setStates(std::vector<unsigned int> states);
		
		unsigned int size() const;
		unsigned int numberOfVariables() const;
		unsigned int getStateIndex(std::vector<unsigned int> _states) const; /*! \brief Returns the index of a state (i.e. 3, 5, 8, 7 => 385992) */

		void setMaxLimit(unsigned int _limit);
		unsigned int getMaxLimit() const;

		void setValues(std::vector<std::vector<unsigned int> > & _values);

	private:
		void calcEnumeration();
		void calcEnumerationValues();

		std::vector<unsigned int> statesPerElement;
		std::vector<std::vector<unsigned int> > enumerations;
		std::vector<std::vector<unsigned int > > values;
		std::vector<std::vector<unsigned int > > enumeratedValues;
		unsigned int combinatorialSize;
		unsigned int maxLimit;
		bool valueSet_flag;
};

inline std::vector<unsigned int> & Enumerator::operator[](size_t _n) {
	return enumerations[_n];
}

inline void Enumerator::setMaxLimit(unsigned int _limit) {
	maxLimit = _limit;
}

inline unsigned int Enumerator::getMaxLimit() const {
	return maxLimit;
}

inline unsigned int Enumerator::size() const {
	return combinatorialSize;
}

inline unsigned int Enumerator::numberOfVariables() const {
	return statesPerElement.size();
}



}

#endif

