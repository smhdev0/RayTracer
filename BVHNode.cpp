//Sam Harry, Student Number 1901522

#include "BVHNode.h"

BVHNode::BVHNode(AABB node, BVHNode * leftNode, BVHNode * rightNode) : node(node), leftNode(leftNode), rightNode(rightNode) {}

BVHNode::BVHNode(AABB node) : node(node), leftNode(NULL), rightNode(NULL) {}

BVHNode::BVHNode() {

}

BVHNode * BVHNode::getNode() {
	return this;
}

AABB BVHNode::getAABBatNode() {
	return this->node;
}

BVHNode * BVHNode::getLeft() {
	return this->leftNode;
}

BVHNode * BVHNode::getRight() {
	return this->rightNode;
}

void BVHNode::setAABBatNode(AABB node) {
	this->node = node;
}

void BVHNode::setLeft(BVHNode * leftNode) {
	this->leftNode = leftNode;
}

void BVHNode::setRight(BVHNode * rightNode) {
	this->rightNode = rightNode;
}

int BVHNode::getDepth() {
	return this->depth;
}

void BVHNode::setDepth(int depth) {
	this->depth = depth;
}

int BVHNode::maxDepth(BVHNode * node) {
	if (node == NULL) {
		return -1;
	}
	else {
		int lDepth = maxDepth(node->getLeft());
		int rDepth = maxDepth(node->getRight());

		if (lDepth > rDepth) {
			return (lDepth + 1);
		}
		else {
			return (rDepth + 1);
		}
	}
}

int BVHNode::findDepth(BVHNode* node, int x) {
	if (node == NULL) {
		return -1;
	}

	int dist = -1;

	if ((node->getAABBatNode().getTriangleList().size() == x)
		|| (dist = findDepth(node->getLeft(), x)) >= 0
		|| (dist = findDepth(node->getRight(), x) >= 0)) {
		
		return dist + 1;

	}
	return dist;
}

BVHNode* BVHNode::getLeftMost(BVHNode * node) {
	while (node->getLeft() != NULL) {
		node = this->getLeft();
	}
	return node;
}