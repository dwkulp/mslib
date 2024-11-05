#ifndef ATOMGRID_H
#define ATOMGRID_H

#include <vector>
#include "AtomPointerVector.h"

/*******************************************************************
   This object creates a 3D grid (based on a given size) and distributes
   atoms in the grid cells

   It can be used to ask for neighbors of an atom i (i.e. the atoms in
   the 27 cells that are neighbors of the atom's cell)

   Useful for algorithms that can use a distance based sub-list of atoms
 *******************************************************************/

namespace MSL { 
class Atom3DGrid {
	public:
		Atom3DGrid();
		Atom3DGrid(AtomPointerVector & _atoms, double _gridSize=1.0);
		~Atom3DGrid();

		bool addAtoms(AtomPointerVector &_ats);
		
		unsigned int getXSize() const;
		unsigned int getYSize() const;
		unsigned int getZSize() const;

		unsigned int size() const;

		AtomPointerVector getCell(unsigned int _i, unsigned int _j, unsigned int _k);
		AtomPointerVector getNeighbors(unsigned int _atomIndex);

		// Get the neighbors for this set of atoms, which are assumed to be in the
		// same coordinate frame.
		bool getNeighbors(AtomPointerVector &_ats,AtomPointerVector &_neighbors, std::string _select="");

		// Check all positions
		void resetAllHiddenFlags();

		// Set a position to skip by adding a 'hidden' flag to each atom with this posId
		void setPositionToSkip(std::string _posId); 
		
	private:
		void setup(AtomPointerVector & _atoms, double _gridSize);
		void buildGrid();
		bool gridBuilt;
		double gridSize;
		AtomPointerVector atoms;
		std::vector<std::vector<unsigned int> > atomIndeces;
		std::vector<std::vector<std::vector<AtomPointerVector> > > grid;
		double xMin;
		double xMax;
		double yMin;
		double yMax;
		double zMin;
		double zMax;
		unsigned int xSize;
		unsigned int ySize;
		unsigned int zSize;

};

inline AtomPointerVector Atom3DGrid::getCell(unsigned int _i, unsigned int _j, unsigned int _k) { return grid[_i][_j][_k]; }
inline unsigned int Atom3DGrid::size() const {return atoms.size();}
inline unsigned int Atom3DGrid::getXSize() const {return xSize;}
inline unsigned int Atom3DGrid::getYSize() const {return ySize;}
inline unsigned int Atom3DGrid::getZSize() const {return zSize;}
}

#endif
