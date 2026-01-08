//Sam Harry, Student Number 1901522

#include <cmath>
#include <iostream>
#include "Vector.h"

Vector::Vector() {

}

Vector::Vector(double i, double j, double k)
	: x(i), y(j), z(k) {}

Vector::~Vector() {

}

	double Vector::magnitude() {
		return sqrt(x * x + y * y + z * z);
	}

	double Vector::getX() {
		return x;
	}

	double Vector::getY() {
		return y;
	}

	double Vector::getZ() {
		return z;
	}

	double Vector::dot(Vector a) {
		return x * a.x + y * a.y + z * a.z;
	}

	Vector Vector::cross(Vector a) {
		return Vector((this->y * a.z - this->z * a.y),
			(this->z * a.x - this->x * a.z),
			(this->x * a.y - this->y * a.x));
	}

	Vector Vector::sub(Vector a) {
		return Vector(x - a.x, y - a.y, z - a.z);
	}

	Vector Vector::add(Vector a) {
		return Vector(x + a.x, y + a.y, z + a.z);
	}

	Vector Vector::mul(double d) {
		return Vector(d * x, d * y, d * z);
	}

	void Vector::normalise() {
		double mag = magnitude();
		if (mag != 0) {
			x = x / mag;
			y = y / mag;
			z = z / mag;
		}
	}

	void Vector::print() {
		std::cout << "x = " << this->x << " , y = " << this->y << " , z = " << this->z << "\n";
	}
