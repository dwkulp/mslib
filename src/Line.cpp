/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2013 The MSL Developer Group (see README.TXT)
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

#include "Line.h"
#include "MslExceptions.h"
#include "PDBFormat.h"

using namespace MSL;
using namespace std;


Line::Line() {

	// Set positional defaults.
	center.setCoor(0.0, 0.0, 0.0);
	direction.setCoor(0.0, 0.0, 0.0);

	// Set color,name,format default
	color.resize(6);
	name         = "line";
	outputFormat = "coord";
}

Line::Line(double _x1, double _y1, double _z1, double _x2, double _y2, double _z2) {

	// Initalize coorinates
	center.setCoor(_x1, _y1, _z1);
	direction.setCoor(_x2, _y2, _z2);

	// Set color,name,format default
	color.resize(6);
	name         = "line";
	outputFormat = "coord";
}
		
Line::Line(const CartesianPoint & _center, const CartesianPoint & _direction) {
	center = _center;
	direction = _direction;

	name         = "line";
	outputFormat = "coord";
}

Line::Line(const Line *_line) {
	center       = _line->center;
	direction    = _line->direction;

	color        = _line->color;
	outputFormat = _line->outputFormat;
	name         = _line->name;
}

Line::Line(const Line& _line) {
	center       = _line.center;
	direction    = _line.direction;
	
	color        = _line.color;
	outputFormat = _line.outputFormat;

}

Line::~Line() {
}


void Line::operator=(const Line& _line) {
	center       = _line.center;
	direction    = _line.direction;
	name         = _line.name;
	color        = _line.color;
	outputFormat = _line.outputFormat;
}

void Line::operator=(Line* _line) {
	center       = _line->center;
	direction    = _line->direction;
	name         = _line->name;
	color        = _line->color;
	outputFormat = _line->outputFormat;
}

bool Line::operator!=(const Line& _line) const {

	if (*this == _line) {
		return false;
	} else {
		return true;
	}
}

bool Line::operator!=(const Line * _line) const {

	if (*this == _line) {
		return false;
	} else {
		return true;
	}

}

bool Line::operator==(const Line &_line) const {

	if (*this == _line) {
		return true;
	} else {
		return false;
	}

	/*
	if (isParallel(_line) && distance(_line) == 0.0) {
		return true;
	} else {
		return false;
	}
	*/
}

bool Line::operator==(const Line *_line) const {

	if (*this == _line) {
		return true;
	} else {
		return false;
	}

	/*
	if (isParallel(*_line) && distance(*_line) == 0.0) {
		return true;
	} else {
		return false;
	}
	*/
}

CartesianPoint Line::getCenter() const {
	return center;
}

void Line::setCenter(const CartesianPoint & _center) {
	center = _center;
}

CartesianPoint Line::getDirection() const {
	return direction;
}

void Line::setDirection(const CartesianPoint & _direction) {
	direction = _direction;
}


double Line::distance(const CartesianPoint & point) const {

	// calculate the distance from a point
	if (direction.length() == 0) {
		throw MslGeneralException("Line::distanceFromPoint, vector has 0 length");
	}

	double d = CartesianGeometry::distanceFromLine(point,center, direction);	
	return d;
}


CartesianPoint Line::projection(const CartesianPoint & point) const {
	// return the point on the line that has minimal distance with a point
	if (direction.length() == 0) {
		throw MslGeneralException("CartesianPoint Line::projection(const CartesianPoint & point) const, vector has 0 length");
	}

	return CartesianGeometry::projection(point,center, center+direction);
	
}

double Line::distance(const Line & _line) const {

	CartesianPoint cx = direction.cross(_line.direction);
	if (cx.length() == 0.0) {
		/********************************************
		 *  case 1, are parallel: all points are equally
		 *  distant from the other line, use the line
		 *  to point distance formula
		 ********************************************/
		return _line.distance(center);
	} else {
		CartesianPoint cd = center - _line.center;
		double d = (cd * cx)/ cx.length();
		if (d<0) {
			return -d;
		} else {
			return d;
		}
	}
	
}

CartesianPoint Line::pointOfMinDistanceFromLine(const Line & _line) const {
	/****************************************************
	 *  point of min distance P' on this line to other line
	 *
	 *             P1
	 *  -----------*---------------------  this line
	 *             |
	 *  -----------*---------------------  _line
	 *             P2
	 ****************************************************/

	vector<CartesianPoint> a = pointsOfMinDistance(_line);
	return a[0];
}

CartesianPoint Line::pointOfMinDistanceToLine(const Line & _line) const {
	/****************************************************
	 *  point of min distance P' on this line to other line
	 *
	 *             P1
	 *  -----------*---------------------  this line
	 *             |
	 *  -----------*---------------------  _line
	 *             P2
	 ****************************************************/
	vector<CartesianPoint> a = pointsOfMinDistance(_line);
	return a[1];
}

vector<CartesianPoint> Line::pointsOfMinDistance(const Line & _line) const {

	/********************************************
	 *  POINT OF MINIMUM DISTANCE BETWEEN TWO LINES
	 *
	 *  We need to find the crossing point of a perpendicular
	 *  line (line3) that intersects line1 and line2.
	 *
	 *  line1 = c1 + t1 * d1
	 *  line2 = c2 + t2 * d2
	 *  line3 = c3 + t3 * d3
	 * 
	 *  where c are the center points, the d are the direction
	 *  vectors of the lines and t are real variables
	 *
	 *  d3 is perpendicular to both and we obtain it as the 
	 *  cross product of d1 and d2
	 *  c3 is an unknown point
	 *
	 *  line1 and line3 cross at p1, line2 and line3 cross at p2,
	 *  we need to find these two points
	 *
	 *  we can set c3 to be the intersecting point on line1:
	 *  c3 = p1 = c1 + t1' * d1    (for a certain t1' value)
	 *
	 *  p2 is a point on line2 and line3 thus
	 *  p2 = c2 + t2' * d2
	 *  p2 = c3 + t3' * d3 = c1 + t1' * d1 + t3' * d3
	 *
	 *  From the above we have that
	 *  c2 + t2' * d2 = c1 + t1' * d1 + t3' * d3
	 *  or
	 *  c1 + t1' * d1 + t3' * d3 - c2 - t2' * d2 = 0
	 *
	 *  Expressed in term of the cartesian coordinates:
	 *  c1 = (xc1, yc1, zc1)
	 *  d1 = (xd1, yd1, zd1)
	 *  c2 = (xc2, yc2, zc2)
	 *  d2 = (xd2, yd2, zd2)
	 *  c3 = (xc3, yc3, zc3)
	 *  d3 = (xd3, yd3, zd3)
	 *
	 *  / xc1 + xd1 * t1' + xd3 * t3' - xc2 - xd2 * t2' = 0
	 * -| yc1 + yd1 * t1' + yd3 * t3' - yc2 - yd2 * t2' = 0
	 *  \ zc1 + zd1 * t1' + zd3 * t3' - zc2 - zd2 * t2' = 0
	 *
	 *  or
	 *
	 *  / xd1 * t1' - xd2 * t2' + xd3 * t3' + xc1 - xc2 = 0
	 * -| yd1 * t1' - yd2 * t2' + yd3 * t3' + yc1 - yc2 = 0
	 *  \ zd1 * t1' - zd2 * t2' + zd3 * t3' + zc1 - zc2 = 0
	 *
	 *  Solving the system for t1', t2', t3', we obtain
	 *  p1 = c1 + t1' * d1
	 *  p2 = c2 + t2' * d2
	 *  
	 *  and we can calculated the p1 p2 distance
	 * 
	 ********************************************/
	CartesianPoint cross = direction.cross(_line.direction);
	vector<CartesianPoint> out;

	if (cross.length() == 0.0) {
		/********************************************
		 *  are parallel, there is an infinite set of
		 *  min distance point, return the center of
		 *  the line and its min distance point on the 
		 *  other line
		 ********************************************/
		 out.push_back(center);
		 out.push_back(_line.projection(center));
		 return out;
	}

	// matrix 3x4
	vector<vector<double> > sys;
	sys.push_back(vector<double>(4, 0.0));
	sys.push_back(vector<double>(4, 0.0));
	sys.push_back(vector<double>(4, 0.0));
	sys[0][0] = direction.getX();
	sys[0][1] = -1 * _line.direction.getX();
	sys[0][2] = cross.getX();
	sys[0][3] = center.getX() - _line.center.getX();
	sys[1][0] = direction.getY();
	sys[1][1] = -1 * _line.direction.getY();
	sys[1][2] = cross.getY();
	sys[1][3] = center.getY() - _line.center.getY();
	sys[2][0] = direction.getZ();
	sys[2][1] = -1 * _line.direction.getZ();
	sys[2][2] = cross.getZ();
	sys[2][3] = center.getZ() - _line.center.getZ();

	vector<double> sol = solve3x3System(sys);
	CartesianPoint p1 = center + direction * sol[0];
	CartesianPoint p2 = _line.center + _line.direction * sol[1];
	//CartesianPoint p3 = p1 + cross * sol[2];

	out.push_back(p1);
	out.push_back(p2);

	return out;
}

double Line::segmentDistance(const Line & _otherSegment) {

	vector<CartesianPoint> v = segmentsClosestPoints(_otherSegment);
	return v[0].distance(v[1]);
}

double Line::segmentDistance(const CartesianPoint & _vec) {
	return CartesianGeometry::distanceFromSegment(_vec, center, direction);
}

double Line::segmentDihedral(const Line & _otherSegment) {
	CartesianPoint p2 = pointOfMinDistanceFromLine(_otherSegment);
	CartesianPoint p1 = p2 + direction;
	CartesianPoint p3 = pointOfMinDistanceToLine(_otherSegment);
	CartesianPoint p4 = p3 + _otherSegment.direction;
	return p1.dihedral(p2, p3, p4);
	
}

vector<CartesianPoint> Line::segmentsClosestPoints(const Line & _otherSegment) {
	/*******************************************************
	 *  The closest points in two segments
	 *   p[0] on this segment
	 *   p[1] on the other segment
	 *******************************************************/
	vector<CartesianPoint> pomd = pointsOfMinDistance(_otherSegment);
	
	// the other end of the segments
	CartesianPoint endOfThis = center + direction;
	CartesianPoint endOfOther = _otherSegment.center + _otherSegment.direction;

	vector<CartesianPoint> out;
	out.push_back(closestPointWithinSegment(pomd[0]));
	out.push_back(_otherSegment.closestPointWithinSegment(pomd[1]));
	return out;
}

CartesianPoint Line::closestPointWithinSegment(CartesianPoint _point) const {

	/*******************************************************
	 *  Given the segment center -> center+direction
	 *  and a point _point, it returns the point on the segment
	 *  that is closest to _point (either _point, if it is within
	 *  the segment, or the closest extremity
	 *******************************************************/
	CartesianPoint p1 = center;
	CartesianPoint p2 = center + direction;
	CartesianPoint pm = projection(_point);
	
	int index = 0;
	for (unsigned int i=0; i<3; i++) {
		if (p1[i] != p2[i]) {
			index = i;
			break;
		}
	}
	
	if (p1[index] < p2[index]) {
		// case 1: the coordinate of the index dimension is greater for p1
		if (pm[index] <= p1[index]) {
			// it is on the p1 side
			return p1;
		} else if (pm[index] >= p2[index]) {
			// it is on the other end side
			return p2;
		} else {
			// it is in between
			return pm;
		}
	} else if (p1[index] > p2[index]) {
		// case 2: the coordinate of the index dimension is greater for p2
		if (pm[index] <= p2[index]) {
			// it is on the p2 side
			return p2;
		} else if (pm[index] >= p1[index]) {
			// it is on the other end side
			return p1;
		} else {
			// it is in between
			return pm;
		}
	} else {
		// segment is a point
		return p1;
	}
}


bool Line::isParallel(const Line & otherLine) const {
	if (direction.cross(otherLine.direction).length() == 0) {
		return true;
	} else {
		return false;
	}
}


bool Line::isOpposite(const Line & otherLine) const {
	if (!isParallel(otherLine) && direction * otherLine.direction < 0.0) {
		return true;
	} else {
		return false;
	}
}


//CartesianPoint Line::pointOfMinDistance(Line point) {
//}

double Line::angle(const CartesianPoint & theVec) const {
	double angle = direction.angle(theVec);
	// make sure it is sharp
	if (angle < -90.0) {
		angle += 180.0;
	}
	if (angle > 90.0) {
		angle -= 180.0;
	}
	return angle;
}

double Line::angle(const Line & otherLine) const {
	// call the angle(CartesianPoint) function with the
	// direction vector of the other line
	return angle(otherLine.direction);
}
		
void Line::translate(const CartesianPoint & vec) {
	center += vec;
}

void Line::xRotate(double degrees, const CartesianPoint & rotationCenter) {

	Matrix m = CartesianGeometry::getXRotationMatrix(degrees);
	center -= rotationCenter;
	center *= m;
	center += rotationCenter;

	direction -= center;
	direction *= m;
	direction += center;
}

void Line::yRotate(double degrees, const CartesianPoint & rotationCenter) {

	Matrix m = CartesianGeometry::getYRotationMatrix(degrees);
	center -= rotationCenter;
	center *= m;
	center += rotationCenter;

	direction -= center;
	direction *= m;
	direction += center;

}

void Line::zRotate(double degrees, const CartesianPoint & rotationCenter) {

	Matrix m = CartesianGeometry::getZRotationMatrix(degrees);
	center -= rotationCenter;
	center *= m;
	center += rotationCenter;

	direction -= center;
	direction *= m;
	direction += center;

}

void Line::rotateAroundAxis(double degrees, const CartesianPoint & rotationCenter, const CartesianPoint & axis) {
	// first rotate the center of the line around the rotation center,
	// then rotate the line direction around the line center

	Matrix m = CartesianGeometry::getRotationMatrix(degrees,axis);
	center -= rotationCenter;
	center *= m;
	center += rotationCenter;

	direction -= center;
	direction *= m;
	direction += center;

}

void Line::rotateAroundAxis(double degrees, const Line & line) {
	// first rotate the center of the line around its projection onto the other line
	// then rotate the line direction around the line center

	CartesianPoint p = projection(center);

	Matrix m = CartesianGeometry::getRotationMatrix(degrees,line.direction);
	center -= p;
	center *= m;
	center += p;

	direction -= center;
	direction *= m;
	direction += center;

}


vector<double> Line::solve3x3System(vector<vector<double> > sys) const {
	vector<double> out(3, 0.0);

	if (sys.size() != 3) {	
		throw MslSizeException("Line::solve3x3System, matrix size invalid variable 'sys', expecting 3");
	}
	if (sys[0].size() != 4) {
		throw MslSizeException("Line::solve3x3System, matrix size invalid variable 'sys[0]',expecting 4");
	}
	if (sys[1].size() != 4) {
		throw MslSizeException("Line::solve3x3System, matrix size invalid variable 'sys[1]',expecting 4");
	}
	if (sys[2].size() != 4) {
		throw MslSizeException("Line::solve3x3System, matrix size invalid variable 'sys[2]',expecting 4");
	}

	//cout << ">>>> " << Matrix(sys).getDeterminant() << endl;

	Matrix m(Matrix(sys).getSubMatrix(0, 2, 0, 2));
	if (m.getDeterminant() == 0) {
		return out;
	}

	/*************************************************************
	 * Swap lines around, so that the biggest number on col [0] is at [0][0] 
	 *************************************************************/
	if (fabs(sys[0][0]) < fabs(sys[1][0]) && fabs(sys[0][0]) < fabs(sys[2][0])) {
		vector<vector<double> >  tmp = sys;
		if (fabs(sys[1][0]) > fabs(sys[2][0])){
			sys[0] = tmp[1];
			sys[1] = tmp[0];
		} else {
			sys[0] = tmp[2];
			sys[2] = tmp[0];
		}
	}

	/*************************************************************
	 * Now normalize the first row so that sys[0][0] = 1
	 *   1  x  x  x
	 *   x  x  x  x
	 *   x  x  x  x
	 *************************************************************/
	double f = sys[0][0];
	for (unsigned int i=0;i<4;i++) {
		sys[0][i] = sys[0][i]/f;
	}


	/*************************************************************
	 * Now zero the first column under sys[0][0]
	 *   1  x  x  x
	 *   0  x  x  x
	 *   0  x  x  x
	 *************************************************************/
	//f = sys[1][0]/sys[0][0];
	f = sys[1][0];
	for (unsigned int i=0;i<4;i++) {
		sys[1][i] = sys[1][i] - sys[0][i] * f;
	}
	//f = sys[2][0]/sys[0][0];
	f = sys[2][0];
	for (unsigned int i=0;i<4;i++) {
		sys[2][i] = sys[2][i] - sys[0][i] * f;
	}

	/*************************************************************
	 * Swap lines around, so that the biggest number on col [0] is at [0][0] 
	 *************************************************************/
	if (fabs(sys[2][1]) > fabs(sys[1][1])) {
		vector<vector<double> >  tmp = sys;
		sys[1] = tmp[2];
		sys[2] = tmp[1];
	}

	
	/*************************************************************
	 * Now normalize the second row so that sys[1][1] = 1
	 *   1  x  x  x
	 *   0  1  x  x
	 *   0  x  x  x
	 *************************************************************/
	f = sys[1][1];
	for (unsigned int i=0;i<4;i++) {
		sys[1][i] = sys[1][i]/f;
	}

	/*************************************************************
	 * Now subtract row 1 * f so that the second column is zeroed
	 *   1  0  x  x
	 *   0  1  x  x
	 *   0  0  x  x
	 *************************************************************/
	f = sys[0][1];
	for (unsigned int i=0;i<4;i++) {
		sys[0][i] = sys[0][i]-sys[1][i] * f;
	}
	f = sys[2][1];
	for (unsigned int i=0;i<4;i++) {
		sys[2][i] = sys[2][i]-sys[1][i] * f;
	}

	/*************************************************************
	 * Now normalize the third row so that aa[2][2] = 1
	 *   1  0  x  x
	 *   0  1  x  x
	 *   0  0  1  x
	 *************************************************************/
	f = sys[2][2];
	for (unsigned int i=0;i<4;i++) {
		sys[2][i] = sys[2][i]/f;
	}

	/*************************************************************
	 * Now subtract row 2 * f so that the third column is zeroed
	 *   1  0  0  x
	 *   0  1  0  x
	 *   0  0  1  x
	 *************************************************************/
	f = sys[0][2];
	for (unsigned int i=0;i<4;i++) {
		sys[0][i] = sys[0][i]-sys[2][i] * f;
	}
	
	f = sys[1][2];
	for (unsigned int i=0;i<4;i++) {
		sys[1][i] = sys[1][i]-sys[2][i] * f;
	}

	/*************************************************************
	 * Once the matrix is diagonalized, finally calculate the solutions
	 * of the system
	 *************************************************************/
	out[0] = -1 * sys[0][3];
	out[1] = -1 * sys[1][3];
	out[2] = -1 * sys[2][3];

	return out;

}

string Line::toString() {

	double axisLength = 1.75;
	double axisWidth  = 0.1;
	if (outputFormat == "pymol"){

		stringstream ss;
		ss << "from pymol.cgo import *"<<endl;
		ss << "from pymol import cmd"<<endl;
		ss << "from pymol.vfont import plain"<<endl;

		if (color.size() < 6){
			color.clear();
			color.push_back(1);
			color.push_back(0);
			color.push_back(0);

			color.push_back(1);
			color.push_back(0);
			color.push_back(0);
		}
		if (name == ""){
			name = "Line";
		}

		// Cylinder +/- axisLength Angstroms from center , in approriate direction
		//		ss << name <<"Direction = [ CYLINDER, "<<center[0] - axisLength*direction[0]<<","<<center[1]- axisLength*direction[1]<<","<<center[2]- axisLength*direction[2]<<","<<(center[0]+axisLength*direction[0])<<","<<(center[1]+axisLength*direction[1])<<","<<(center[2]+axisLength*direction[2])<< ","<<axisWidth<<", "<<color[0]<<","<<color[1]<<","<<color[2]<<","<<color[3]<<","<<color[4]<<","<<color[5]<<" ]"<<endl;
		ss << name <<"Direction = [ CYLINDER, "<<center[0]<<","<<center[1]<<","<<center[2]<<","<<(center[0]+axisLength*direction[0])<<","<<(center[1]+axisLength*direction[1])<<","<<(center[2]+axisLength*direction[2])<< ","<<axisWidth<<", "<<color[0]<<","<<color[1]<<","<<color[2]<<","<<color[3]<<","<<color[4]<<","<<color[5]<<" ]"<<endl;
		ss << "cmd.load_cgo("<<name<<"Direction,'"<<name<<"Direction');"<<endl;


		return ss.str();

	} else if (outputFormat == "pymolAtom") {

		stringstream ss;
		ss << "from pymol.cgo import *"<<endl;
		ss << "from pymol import cmd"<<endl;
		ss << "from pymol.vfont import plain"<<endl;

		if (color.size() < 6){
			color.clear();
			color.push_back(1);
			color.push_back(0);
			color.push_back(0);

			color.push_back(1);
			color.push_back(0);
			color.push_back(0);
		}
		if (name == ""){
			name = "Line";
		}
		// Cylinder +/- axisLength Angstroms from center , in approriate direction
		ss << name <<"Direction = [ CYLINDER, "<<center[0] - axisLength*direction[0]<<","<<center[1]- axisLength*direction[1]<<","<<center[2]- axisLength*direction[2]<<","<<(center[0]+axisLength*direction[0])<<","<<(center[1]+axisLength*direction[1])<<","<<(center[2]+axisLength*direction[2])<< ","<<axisWidth<<", "<<color[0]<<","<<color[1]<<","<<color[2]<<","<<color[3]<<","<<color[4]<<","<<color[5]<<" ]"<<endl;
		ss << "cmd.load_cgo("<<name<<"Direction,'"<<name<<"Direction');"<<endl;

		
		Real x = center[0] + axisLength*direction[0];
		Real y = center[1] + axisLength*direction[1];
		Real z = center[2] + axisLength*direction[2];
		
		string pdbline = PDBFormat::createAtomLine(PDBFormat::createAtomData((string)"DUM", x,y,z,(string)"H"));

		string pdbName = name+"Atom";

		ss << "cmd.read_pdbstr(\""<<pdbline<<"\",\""<<pdbName<<"\")"<<endl;
		
		return ss.str();
	} else {
		char c [200];
		sprintf(c, "[%10.3f %10.3f %10.3f],[%10.3f %10.3f %10.3f] ", center[0],center[1],center[2],direction[0],direction[1],direction[2]);
		return (string)c;
	}

}

void Line::setColor(double _red, double _green, double _blue){
  if (color.size() < 6) {
	  color.resize(6);
  } 
  color[0] = _red; 
  color[1] = _blue; 
  color[2] = _green;
  color[3] = _red; 
  color[4] = _blue; 
  color[5] = _green;
}

void Line::setColor(double _red, double _green, double _blue, double _red2, double _green2, double _blue2){
  if (color.size() < 6) {
	  color.resize(6);
  } 
  color[0] = _red; 
  color[1] = _blue; 
  color[2] = _green;
  color[3] = _red2; 
  color[4] = _blue2; 
  color[5] = _green2;
}


