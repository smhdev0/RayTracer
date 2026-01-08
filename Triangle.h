//Sam Harry, Student Number 1901522

#pragma once
#include "Point.h"
#include "Triangle.h"

class AABB;

class Triangle {
	Point a, b, c;

public:

	Point getA();

	Point getB();
	
	Point getC();

	Triangle();

	Triangle(Point a, Point b, Point c);

	~Triangle();

	void print();

	double getMaxXCoord();

	double getMinXCoord();

	double getMaxYCoord();

	double getMinYCoord();

	double getMaxZCoord();

	double getMinZCoord();

	double getCentreX();

	double getCentreY();

	double area(double x1, double y1, double x2, double y2, double x3, double y3);

	bool isInBounds(AABB * box);

};