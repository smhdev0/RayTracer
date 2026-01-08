//Sam Harry, Student Number 1901522

#pragma once

class Vector {
	double x, y, z;

public :

	Vector();

	Vector(double x, double y, double z);

	~Vector();

	double magnitude();

	double getX();

	double getY();

	double getZ();

	double dot(Vector);

	Vector cross(Vector);

	Vector sub(Vector);

	Vector add(Vector);

	Vector mul(double);

	void normalise();

	void print();
};