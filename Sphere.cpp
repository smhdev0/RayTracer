//Sam Harry, Student Number 1901522

#include "Point.h"
#include "Sphere.h"
#include <vector>

Sphere::Sphere() {

}

Sphere::Sphere(double radius, Vector sphereCentre, double rCol, double gCol, double bCol)
	: r(radius), cs(sphereCentre), colR(rCol), colG(gCol), colB(bCol) {}


Sphere::~Sphere() {

}

void Sphere::intersect() const {
}
