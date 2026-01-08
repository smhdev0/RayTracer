//Sam Harry, Student Number 1901522

#pragma once
#include "Triangle.h"
#include <vector>

class AABB {

	Point bl;
	Point tr;
	std::vector<Triangle> triangleList;
	bool hasTested;
	bool isSplit = false;

public :
	AABB();

	AABB(Point bl, Point tr);

	~AABB();

	void intersect();

	double getLowestXBound();

	double getHighestYBound();

	double getHighestXBound();

	double getLowestYBound();

	double getWidth();

	double getHeight();

	Point getBottomLeft();

	Point getTopRight();

	void setBottomLeft(Point bl);

	void setTopRight(Point tr);

	void setTriangleList(std::vector<Triangle> triangleList);

	void setTested();

	void setSplit();

	bool getSplit();

	bool getTested();

	std::vector<Triangle> & getTriangleList();
};