//Sam Harry, Student Number 1901522

#pragma once

class Point{
	double x, y, z;

public:

	Point();

	Point(double x, double y, double z);

	double getX();

	double getY();

	double getZ();

	~Point();

	Point sub(Point a);

	void print();
};