//Sam Harry, Student Number 1901522

#include "Point.h"
#include "Triangle.h"
#include "AABB.h"
#include <iostream>
#include <algorithm>
#include "Vector.h"

Triangle::Triangle() {

}

Triangle::Triangle(Point x, Point y, Point z) 
	: a(x), b(y), c(z) {}


Triangle::~Triangle() {

}

Point Triangle::getA() {
	return this->a;
}

Point Triangle::getB() {
	return this->b;
}

Point Triangle::getC() {
	return this->c;
}

void Triangle::print() {
	a.print();
	b.print();
	c.print();
}

double Triangle::getMaxXCoord() {
	return std::max({ a.getX(), b.getX(), c.getX() });
}

double Triangle::getMinXCoord() {
	return std::min({ a.getX(), b.getX(), c.getX() });
}

double Triangle::getMaxYCoord() {
	return std::max({ a.getY(), b.getY(), c.getY() });
}

double Triangle::getMinYCoord() {
	return std::min({ a.getY(), b.getY(), c.getY() });
}

double Triangle::getMaxZCoord() {
	return std::max({ a.getZ(), b.getZ(), c.getZ() });
}

double Triangle::getMinZCoord() {
	return std::min({ a.getZ(), b.getZ(), c.getZ() });
}

double Triangle::getCentreX() {
	return (a.getX() + b.getX() + c.getX()) / 3;
}

double Triangle::getCentreY() {
	return (a.getY() + b.getY() + c.getY()) / 3;
}

double Triangle::area(double x1, double y1, double x2, double y2, double x3, double y3) {
	return abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2);
}

bool Triangle::isInBounds(AABB* box) {

	if ((this->getMinXCoord() <= box->getHighestXBound() && this->getMaxXCoord() >= box->getLowestXBound()) &&
		(this->getMinYCoord() <= box->getHighestYBound() && this->getMaxYCoord() >= box->getLowestYBound())) {
		return true;
	}

	return false;
}

