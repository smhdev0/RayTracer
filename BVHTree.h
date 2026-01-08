//Sam Harry, Student Number 1901522

#pragma once
#include "BVHNode.h"
#include "AABB.h"
class BVHTree
{
	BVHNode * root;

public :
	BVHTree(BVHNode * root);

	void setRoot(BVHNode * root);

	BVHNode * getRoot();


};

