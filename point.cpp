//Sam Harry, Student Number 1901522

#include <cmath>
#include <iostream>
#include "Point.h"

Point::Point() {

}

Point::Point(double i, double j, double k)
	: x(i), y(j), z(k) {}

double Point::getX() {
	return this->x;
}

double Point::getY() {
	return this->y;
}

double Point::getZ() {
	return this->z;
}

Point::~Point() {

}

Point Point::sub(Point a) {
	return Point(x - a.x, y - a.y, z - a.z);
}


void Point::print() {
	std::cout << "x = " << this->x << " y = " << this->y << " x = " << this->z << std::endl;
}

	
