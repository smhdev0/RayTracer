//Sam Harry, Student Number 1901522

#pragma once
#include "Vector.h"
#include <vector>

class Sphere {
	

public :

	Vector cs;

	double r;

	double colR;

	double colG;
	
	double colB;

	Sphere();

	Sphere(double r, Vector cs, double colR, double colG, double colB);

	void intersect() const;

	~Sphere();
};