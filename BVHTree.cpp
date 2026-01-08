//Sam Harry, Student Number 1901522

#include "BVHTree.h"

BVHTree::BVHTree(BVHNode * root) : root(root) {}

void BVHTree::setRoot(BVHNode * root) {
	this->root = root;
}

BVHNode * BVHTree::getRoot() {
	return this->root;
}
