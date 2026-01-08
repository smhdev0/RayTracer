//Sam Harry, Student Number 1901522

#include "AABB.h"
#include <vector>
#include "Triangle.h"
#include "Point.h"

AABB::AABB() {}

AABB::AABB(Point bl, Point tr) : bl(bl), tr(tr) {}

AABB::~AABB() {}

void AABB::intersect() {}

double AABB::getLowestXBound() {
	return this->bl.getX();
}

double AABB::getHighestYBound() {
	return this->tr.getY();
}

double AABB::getHighestXBound() {
	return this->tr.getX();
}

double AABB::getLowestYBound() {
	return this->bl.getY();
}

double AABB::getWidth() {
	return this->getHighestXBound() - this->getLowestXBound();
}

double AABB::getHeight() {
	return this->getHighestYBound() - this->getLowestYBound();
}

Point AABB::getBottomLeft() {
	return this->bl;
}

void AABB::setTested() {
	this->hasTested = true;
}

void AABB::setBottomLeft(Point bl) {
	this->bl = bl;
}

void AABB::setTopRight(Point tr) {
	this->tr = tr;
}

bool AABB::getTested() {
	return this->hasTested;
}

void AABB::setSplit() {
	this->isSplit = true;
}

bool AABB::getSplit() {
	return this->isSplit;
}


Point AABB::getTopRight() {
	return this->tr;
}

void AABB::setTriangleList(std::vector<Triangle> triangleList) {
	this->triangleList = triangleList;
}

std::vector<Triangle> & AABB::getTriangleList() {
	return this->triangleList;
}

