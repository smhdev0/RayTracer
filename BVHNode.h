//Sam Harry, Student Number 1901522

#pragma once
#include "AABB.h"
class BVHNode
{
	AABB node;
	BVHNode * leftNode;
	BVHNode* rightNode;
	int depth;

public:

	BVHNode(AABB node, BVHNode * leftNode, BVHNode * rightNode);

	BVHNode(AABB node);

	BVHNode();
	
	void setAABBatNode(AABB node);

	void setLeft(BVHNode * node);

	void setRight(BVHNode * node);

	AABB getAABBatNode();

	BVHNode * getLeft();

	BVHNode * getRight();

	BVHNode * getNode();

	int maxDepth(BVHNode* node);

	int findDepth(BVHNode* node, int x);

	int getDepth();

	void setDepth(int depth);

	BVHNode* getLeftMost(BVHNode * node);
	
};

