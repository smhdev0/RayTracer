//Sam Harry, Student Number 1901522

#include "cMain.h"
#include "Vector.h"
#include "point.h"
#include "plyRead.h"
#include "Triangle.h"
#include "plyRead.h"
#include "Sphere.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <list>
#include <thread>
#include <omp.h>
#include "AABB.h"
#include "BVHTree.h"

wxBEGIN_EVENT_TABLE(cMain, wxFrame)
EVT_BUTTON(10001, OnButtonClicked)
EVT_PAINT(OnPaint)

wxEND_EVENT_TABLE()

int boxSize = 0;
int boxRowCols = 4;
int hRes = 1000; //horizontal resolution of the image
int vRes = 1000; //vertical resolution of the image
std::vector<AABB> leafBoxList; //global vector containing the list of leaves to trace.

//Method to ray trace a triangle mesh without an acceleration structure
void cMain::RayTraceTriangleMesh(wxImage* graphic, std::vector<Triangle> & triangleList) {

#pragma omp parallel
	{

		Triangle triangle; //current triangle to test
		Vector reflection;
		Vector origin;
		Vector light = Vector(20, 50, 30);
		Vector Lnorm;
		Vector d = Vector(0, 0, 1);
		double t; //intersection point
		double D; //distance from origin to the plane
		Vector pHit; //point on the plane that has been hit
		Vector planeNorm; //plane normal for the triangle
		double rCol;
		double gCol;
		double bCol;
		int imageWidth = hRes;
		int imageHeight = vRes;
		double colDiffuse;
		double colPhong;
		double phong = 64;
		double foreTriangleZ; //triangle closest to the camera
		double imageAspectRatio = imageWidth / vRes;
		int fovDegrees = 25;
		double fovRadians = fovDegrees * (atan(1) * 4) / 180; //convert fov to radians

		#pragma parallel for collapse(3)
		for (int i = 0; i < imageWidth; i++) {
			for (int j = 0; j < imageHeight; j++) {
				foreTriangleZ = 1;
				for (int k = 0; k < triangleList.size(); k++) {
					triangle = triangleList.at(k);

					//update the origin for every pixel

					 //direction, +ve x means away from the camera, -ve towards the camera.
					double x = .60 * (2 * ((i + 0.5) / imageWidth) - 1) * imageAspectRatio * tan(fovRadians / 2);
					double y = .6 * (((1 - 2 * ((j + 0.5) / imageHeight)) * tan(fovRadians / 2)) + 0.17);

					origin = Vector(x, y, -1);
					//colPhong = 0.8;

					origin.normalise();

					rCol = 1;
					gCol = 1;
					bCol = 1;

					//Step 1 - finding the ray hit point

					//ray equation --------> P = O +tR
					//P = intersection point, O = origin, t scalar value of length along the ray, R = direction of the ray

					//plane equation ----------> Ax + By + Cz + D = 0
					//A, B, C are the coordinates of the normal to the plane
					//D is the distance from the origin to the plane.
					Vector a = Vector(triangle.getA().getX(), triangle.getA().getY(), triangle.getA().getZ());

					Vector b = Vector(triangle.getB().getX(), triangle.getB().getY(), triangle.getB().getZ());

					Vector c = Vector(triangle.getC().getX(), triangle.getC().getY(), triangle.getC().getZ());

					Vector ba = a.sub(b);

					Vector cb = b.sub(c);

					Vector ca = a.sub(c);

					planeNorm = ba.cross(cb); //plane normal is the cross product of two sides of the triangle


					planeNorm.normalise();

					//is the ray parallel to the triangle? Test if the plane normal is perpendicular to the ray direction
					if (planeNorm.dot(d) != 0) {

						//backface culling
						if (planeNorm.dot(d) < 0) {
							continue;
						}

						//compute D by subbing in one of the triangle vertices to the plane equation
						D = planeNorm.dot(Vector(triangle.getA().getX(), triangle.getA().getY(), triangle.getA().getZ()));
						//substitute P from the ray equation into the plane equation:
						t = -(planeNorm.dot(origin) - D) / (planeNorm.dot(d));

						if (t < 0) {
							continue;
						}

						//calculate the position of P using the ray equation:
						pHit = origin.add(d.mul(t));

						//Step 2 - is the ray inside or outside the triangle?
						Vector C;

						//test ab
						Vector vp0 = pHit.sub(b);
						C = ba.cross(vp0);

						if (planeNorm.dot(C) > 0) {
							continue;
						}

						Vector vp1 = pHit.sub(c);
						C = cb.cross(vp1);

						if (planeNorm.dot(C) > 0) {
							continue;
						}

						Vector vp2 = pHit.sub(c);
						C = ca.cross(vp2);

						if (planeNorm.dot(C) < 0) {
							continue;
						}


						Lnorm = light.sub(pHit); //light normal calculated by light vector subtract intersection vector

						Lnorm.normalise();

						colDiffuse = planeNorm.dot(Lnorm);

						if (colDiffuse <= 0) colDiffuse = 0; //if costheta is less that 0, it is in shadow
						if (colDiffuse > 1) colDiffuse = 1; //if not, use the costheta value

						rCol = .1 + colDiffuse;//total colour in each channel is the ambient + diffuse + phong, all multiplied by the light colour
						gCol = .1 + colDiffuse;
						bCol = .1 + colDiffuse;

						if (rCol > 1) rCol = 1;
						if (rCol < 0) rCol = 0;
						if (gCol > 1) gCol = 1;
						if (gCol < 0) gCol = 0;
						if (bCol > 1) bCol = 1;
						if (bCol < 0) bCol = 0;

						//rCol = 1;
						//gCol = 1;
						//bCol = 1;
						//Step 3 
					}
					graphic->SetRGB(i, j, rCol * 255, gCol * 255, bCol * 255);

				}

			}
		}
	}
}

//Method to ray trace a triangle mesh using an acceleration structure (Bounding Volume Hierarchy)
void cMain::RayTraceBVH(wxImage* graphic, std::vector<Triangle>& triangleList, std::vector<AABB>& AABBlist, BVHTree & tree) {

	omp_set_num_threads(6);
#pragma omp parallel
	{

		Triangle triangle; //current triangle to test

		Vector reflection;
		Vector origin;
		Vector light = Vector(20, 50, 30);
		Vector Lnorm;
		Vector d = Vector(0, 0, 1);
		double t; //intersection point
		double D; //distance from origin to the plane
		Vector pHit; //point on the plane that has been hit
		Vector planeNorm; //plane normal for the triangle
		double rCol;
		double gCol;
		double bCol;
		int imageWidth = hRes;
		int imageHeight = vRes;
		double colDiffuse;
		double colPhong;
		double phong = 64;
		double foreTriangleZ; //triangle closest to the camera
		double imageAspectRatio = imageWidth / vRes;
		int fovDegrees = 25;
		double fovRadians = fovDegrees * (atan(1) * 4) / 180; //convert fov to radians
		bool hasHit; 
		AABB testBox;
		AABB baseBox = AABBintersect(triangleList);

		#pragma parallel for collapse(3)
		for (int i = 0; i < imageWidth; i++) {
			for (int j = 0; j < imageHeight; j++) {
				hasHit = false;
				foreTriangleZ = 1;
				double x = .60 * (2 * ((i + 0.5) / imageWidth) - 1) * imageAspectRatio * tan(fovRadians / 2);
				double y = .60 * (((1 - 2 * ((j + 0.5) / imageHeight)) * tan(fovRadians / 2)) + 0.17);


				//firstly, write off any rays that do not intersect the base bounding box.
				if (x < baseBox.getLowestXBound()) {
					continue;
				}

				if (y < baseBox.getLowestYBound()) {
					continue;

				}

				if (x > baseBox.getHighestXBound()) {
					continue;

				}

				if (y > baseBox.getHighestYBound()) {
					continue;
				}

				for (int b = 0; b < leafBoxList.size(); b++) {
					testBox = leafBoxList.at(b);

					//Check if the calculated world coordinate for the ray would enter the bounding box

					if (x < testBox.getLowestXBound()) {
						continue;
					}

					if (y < testBox.getLowestYBound()) {
						continue;

					}

					if (x > testBox.getHighestXBound()) {
						continue;

					}

					if (y > testBox.getHighestYBound()) {
						continue;
					}

					else {
						hasHit = true;
						break;
					}
				}

				if (!hasHit) {
					continue;
				}

				
				
			//CODE FOR TRAVERSING THE TREE - UNUSED AS INEFFICIENT
			//testBox = getAABBFromTree(tree, x, y);
			//if (&testBox == nullptr) {
			//	continue;
			//}

				//BVHNode* hitNode = tree.getRoot();
				
				//BVHNode* hitNode = getNodeHit(tree.getRoot(), x, y);

				//if (hitNode->getAABBatNode().getTriangleList().size() == baseBox.getTriangleList().size()) {
				//	continue;
				//}
				
				//testBox = hitNode->getAABBatNode();

				//if (testBox.getHighestXBound() == x || testBox.getLowestXBound() == x + 1) {
				//	graphic->SetRGB(i, j, 255, 255, 255);
				//}

				//if (testBox.getHighestYBound() == y || testBox.getLowestYBound() == y) {
				//	graphic->SetRGB(i, j, 255, 255, 255);
				//}


					for (int k = 0; k < testBox.getTriangleList().size(); k++) {
						triangle = testBox.getTriangleList().at(k);

						//update the origin for every pixel

						origin = Vector(x, y, -1);

						rCol = 0;
						gCol = 0;
						bCol = 0;

						//Step 1 - finding the ray hit point

						//ray equation --------> P = O +tR
						//P = intersection point, O = origin, t scalar value of length along the ray, R = direction of the ray

						//plane equation ----------> Ax + By + Cz + D = 0
						//A, B, C are the coordinates of the normal to the plane
						//D is the distance from the origin to the plane.
						Vector a = Vector(triangle.getA().getX(), triangle.getA().getY(), triangle.getA().getZ());

						Vector b = Vector(triangle.getB().getX(), triangle.getB().getY(), triangle.getB().getZ());

						Vector c = Vector(triangle.getC().getX(), triangle.getC().getY(), triangle.getC().getZ());

						Vector ba = a.sub(b);

						Vector cb = b.sub(c);

						Vector ca = a.sub(c);

						planeNorm = ba.cross(cb); //plane normal is the cross product of two sides of the triangle

						planeNorm.normalise();

						//is the ray parallel to the triangle? Test if the plane normal is perpendicular to the ray direction
						if (planeNorm.dot(d) != 0) {

							//backface culling
							if (planeNorm.dot(d) < 0) {
								continue;
							}

							//compute D by subbing in one of the triangle vertices to the plane equation
							D = planeNorm.dot(Vector(triangle.getA().getX(), triangle.getA().getY(), triangle.getA().getZ()));
							//substitute P from the ray equation into the plane equation:
							t = -(planeNorm.dot(origin) - D) / (planeNorm.dot(d));

							if (t < 0) {
								continue;
							}

							//calculate the position of P using the ray equation:
							pHit = origin.add(d.mul(t));

							//Step 2 - is the ray inside or outside the triangle?
							Vector C;

							//test ab
							Vector vp0 = pHit.sub(b);
							C = ba.cross(vp0);

							if (planeNorm.dot(C) > 0) {
								continue;
							}

							Vector vp1 = pHit.sub(c);
							C = cb.cross(vp1);

							if (planeNorm.dot(C) > 0) {
								continue;
							}

							Vector vp2 = pHit.sub(c);
							C = ca.cross(vp2);

							if (planeNorm.dot(C) < 0) {
								continue;
							}

							Lnorm = light.sub(pHit); //light normal calculated by light vector subtract intersection vector

							Lnorm.normalise();

							colDiffuse = planeNorm.dot(Lnorm);

							if (colDiffuse <= 0) colDiffuse = 0; //if costheta is less that 0, it is in shadow
							if (colDiffuse > 1) colDiffuse = 1; //if not, use the costheta value


							rCol = .1 + colDiffuse; //total colour in each channel is the ambient + diffuse + phong, all multiplied by the light colour
							gCol = .1 + colDiffuse;
							bCol = .1 + colDiffuse;

							if (rCol > 1) rCol = 1;
							if (rCol < 0) rCol = 1;
							if (gCol > 1) gCol = 1;
							if (gCol < 0) gCol = 1;
							if (bCol > 1) bCol = 1;
							if (bCol < 0) bCol = 1;

							//Step 3 
						}
						//set the graphic to the established colour
						
						graphic->SetRGB(i, j, rCol * 255, gCol * 255, bCol * 255);
						
					}

			}
			//print out the current row we are rendering (shows the progression of rendering)
			//std::cout << "i " << i << std::endl;
		}
		}

	}
	
	//funtion to check whether a ray is within a bounding volume's x and y axis.
	bool cMain::isWithinAABB(AABB * aabb, double x, double y) {

		if (x < aabb->getLowestXBound()) {
			return false;
		}

		if (y < aabb->getLowestYBound()) {
			return false;

		}

		if (x > aabb->getHighestXBound()) {
			return false;
		}

		if (y > aabb->getHighestYBound()) {
			return false;
		}

		return true;
	}

	//used to populate the leaf box vector, as tree traversal is inefficient
	void cMain::populateBoxList(BVHNode* node) {

		if (node->getLeft() != NULL) {
			populateBoxList(node->getLeft());
		}

		if (node->getRight() != NULL) {
			populateBoxList(node->getRight());
		}

		if (node->getLeft() == NULL && node->getRight() == NULL) {
			leafBoxList.push_back(node->getAABBatNode());
		}
	}

	//used to traverse the hierarchy and reach the leaf that has been hit by a ray (inefficient)
	BVHNode* cMain::getNodeHit(BVHNode* node, double x, double y) {

			if (node->getLeft() != NULL) {
				AABB box = node->getAABBatNode();
				if (x > box.getLowestXBound() &&
					y > box.getLowestYBound() &&
					x < box.getHighestXBound() &&
					y < box.getHighestYBound()) {
					getNodeHit(node->getLeft(), x, y);
				}
				
			}

			if (node->getRight() != NULL) {
				AABB box = node->getAABBatNode();
				if (x > box.getLowestXBound() &&
					y > box.getLowestYBound() &&
					x < box.getHighestXBound() &&
					y < box.getHighestYBound()) {
					getNodeHit(node->getRight(), x, y);
				}

			}
				return node;

	}


	//A method that I created when experimenting with splitting the image into boxes before the BVH implementation
	std::vector<AABB> cMain::testSplit(std::vector<Triangle>& triangleList) {
		std::vector<AABB> AABBlist;
		AABB baseBox = AABBintersect(triangleList);

		std::cout << "Top right : " << baseBox.getTopRight().getX();
		
		double middleX = baseBox.getTopRight().getX() / 2;

		std::cout << "Middel x: " << middleX << std::endl;
		Point middleTop = Point(middleX, baseBox.getTopRight().getY(), baseBox.getTopRight().getZ());

		middleTop.print();

		Point middleBottom = Point(middleX, baseBox.getBottomLeft().getY(), baseBox.getTopRight().getZ());

		middleBottom.print();

		AABB leftBox = AABB(baseBox.getBottomLeft(), middleTop);
		AABB rightBox = AABB(middleBottom, baseBox.getTopRight());
		std::vector<Triangle> newTriangleListLeft;
		std::vector<Triangle> newTriangleListRight;

		for (int k = 0; k < triangleList.size(); k++) {
			Triangle triangle = triangleList.at(k);

			if (triangle.isInBounds(&leftBox)) {
				newTriangleListLeft.push_back(triangle);
			}

			if (triangle.isInBounds(&rightBox)) {
				newTriangleListRight.push_back(triangle);

			}
		}

		leftBox.setTriangleList(newTriangleListLeft);
		rightBox.setTriangleList(newTriangleListRight);

		std::cout << newTriangleListLeft.size() << std::endl;
		std::cout << newTriangleListRight.size() << std::endl;


		AABBlist.push_back(leftBox);
		AABBlist.push_back(rightBox);

		return AABBlist;
		
	}

	//A method to create a list of triangles for a given bounding box, for use when creating the 
	//bounding volume hierarchy
	std::vector<Triangle> cMain::getAABBTriangles(AABB * box, std::vector<Triangle> & triangleList) {

		std::vector <Triangle> newTriangleList;
		for (int i = 0; i < triangleList.size(); i++) {
			if (triangleList.at(i).isInBounds(box)) {
				newTriangleList.push_back(triangleList.at(i));
			}
		}

		return newTriangleList;
		
	}

	//get the average middle of the triangles for splitting
	double cMain::getAverageTriangleCentroid(AABB* box, bool isVerticalSplit) {
		std::vector <Triangle> boxTriangleList = box->getTriangleList();

		if (isVerticalSplit) {
			double sumX = 0;
			for (int i = 0; i < boxTriangleList.size(); i++) {
				sumX += boxTriangleList.at(i).getCentreX();
			}

			return sumX / boxTriangleList.size();
		}

		double sumY = 0;

		for (int i = 0; i < boxTriangleList.size(); i++) {
			sumY += boxTriangleList.at(i).getCentreY();
		}

		return sumY / boxTriangleList.size();
		
	}

	//function to recursively construct a bounding volume hierarchy
	void cMain::recursiveConstructBVHBox(BVHNode * node, BVHTree * tree, bool isVerticalSplit) {

		if (node == NULL) {
			return;
		}

		//uncomment to use a recursive depth
		//if (node->getDepth() > 200) {
		//	return;
		//}

		AABB box = node->getNode()->getAABBatNode();

		double middleX = box.getTopRight().getX() / 2;
		double middleY = box.getTopRight().getY() / 2;
		AABB boxOne;
		AABB boxTwo;

		
		if (isVerticalSplit) {
			double splitPointX = getAverageTriangleCentroid(&box, isVerticalSplit);
			Point middleTop = Point(splitPointX, box.getTopRight().getY(), box.getTopRight().getZ());
			Point middleBottom = Point(splitPointX, box.getBottomLeft().getY(), box.getTopRight().getZ());
			boxOne = AABB(box.getBottomLeft(), middleTop);
			boxTwo = AABB(middleBottom, box.getTopRight());
		}
		else {
			double splitPointY = getAverageTriangleCentroid(&box, isVerticalSplit);
			Point middleRight = Point(box.getTopRight().getX(), splitPointY, box.getTopRight().getZ());
			Point middleLeft = Point(box.getBottomLeft().getX(), splitPointY, box.getTopRight().getZ());
			boxOne = AABB(box.getBottomLeft(), middleRight);
			boxTwo = AABB(middleLeft, box.getTopRight());
		}

		boxOne.setTriangleList(getAABBTriangles(&boxOne, box.getTriangleList()));
		//std::cout << box.getTriangleList().size() << std::endl;
		boxTwo.setTriangleList(getAABBTriangles(&boxTwo, box.getTriangleList()));
		//std::cout << boxTwo.getTriangleList().size() << std::endl;

		
		if (boxOne.getTriangleList().size() != 0) {
			BVHNode* leftNode = new BVHNode(boxOne);
			leftNode->setDepth(node->getDepth() + 1);
			node->setLeft(leftNode);
			boxSize++;

		}

		if (boxTwo.getTriangleList().size() != 0) {
			BVHNode* rightNode = new BVHNode(boxTwo);
			rightNode->setDepth(node->getDepth() + 1);
			node->setRight(rightNode);
			boxSize++;

		}

		//if the box has less than 200 triangles, return from the recursion
		if (boxOne.getTriangleList().size() > 200) {
			recursiveConstructBVHBox(node->getLeft(), tree, !isVerticalSplit);

		}
		
		if (boxTwo.getTriangleList().size() > 200) {
			recursiveConstructBVHBox(node->getRight(), tree, !isVerticalSplit);
		}
		return;

	}

	//call to the recursive function
	void cMain::recursiveConstructBVH(BVHTree * tree, bool isVerticalSplit) {

		recursiveConstructBVHBox(tree->getRoot(), tree, isVerticalSplit);
		std::cout << "Number of boxes: " << boxSize << std::endl;

	}
				

	//A method I created to test how the bounding volume hierarchy would work. This is not a 
	//hierarchy, but rather a grid with no recursion.
	std::vector<AABB> cMain::constructBVH(std::vector<Triangle>& triangleList, wxImage & graphic) {

		AABB baseBox = AABBintersect(triangleList);

		std::vector<AABB> AABBlist;

				for (int i = 0; i < 20; i++) {
					for (int j = 0; j < 20; j++) {
						std::vector<Triangle> newTriangleList;

						Point bottomLeft = Point((baseBox.getLowestXBound() + i * (baseBox.getWidth() / 20)), (baseBox.getLowestYBound() + j * (baseBox.getHeight() / 20)), baseBox.getBottomLeft().getZ());
						Point topRight = Point((bottomLeft.getX() + (baseBox.getWidth() / 20)), (bottomLeft.getY() + (baseBox.getHeight() / 20)), baseBox.getBottomLeft().getZ());

						AABB newAABB = AABB(bottomLeft, topRight);


						for (int k = 0; k < triangleList.size(); k++) {
							Triangle triangle = triangleList.at(k);

 							if (triangle.isInBounds(&newAABB)) {
								newTriangleList.push_back(triangle);
							}

						}

						if (newTriangleList.size() == 0) {
							continue;
						} 

						

						newAABB.setTriangleList(newTriangleList);

						if (newAABB.getTriangleList().size() != 0) {
							AABBlist.push_back(newAABB);
						}
					}
				}
				
	return AABBlist;
	
	}


	//A method used to create a base bounding box for a given triangle list
	AABB cMain::AABBintersect(std::vector<Triangle> &triangleList) {

		double highestX = 0, lowestX = 0, highestY = 0, lowestY = 0, highestZ = -1, lowestZ = 1;

		for (int k = 0; k < triangleList.size(); k++) {
			Triangle triangle = triangleList.at(k);

			if (triangle.getMaxXCoord() > highestX) {
				highestX = triangle.getMaxXCoord();
			}

			if (triangle.getMinXCoord() < lowestX) {
				lowestX = triangle.getMinXCoord();
			}

			if (triangle.getMaxYCoord() > highestY) {
				highestY = triangle.getMaxYCoord();
			}

			if (triangle.getMinYCoord() < lowestY) {
				lowestY = triangle.getMinYCoord();
			}

			if (triangle.getMaxZCoord() > highestZ) {
				highestZ = triangle.getMaxZCoord();
			}

			if (triangle.getMinZCoord() < lowestZ) {
				lowestZ = triangle.getMinZCoord();
			}
		}

		Point bottomLeft = Point(lowestX, lowestY, lowestZ);
		Point topRight = Point(highestX, highestY, lowestZ);

		AABB box = AABB(bottomLeft, topRight);
		box.setTriangleList(triangleList);

		return box; 
	}


//A method to render spheres with diffuse and phong reflection
void cMain::RayTraceSphere(std::vector<Sphere> *sphereList, wxImage &graphic) {

	//3d line = p = o + dt
	//sphere = (p - c)squared = r squared
	//find where these two intersect

	/**
	int o, //origin of ray (3d coordinate)
		d, //direction of ray
		c, //centre of sphere
		r, //radius of sphere
		p, //3d points on sphere
		t; //solution we solve
	**/

	int imageWidth = hRes;
	int imageHeight = vRes;
	double imageAspectRatio = imageWidth / imageHeight;
	int fovDegrees = 90;
	double fovRadians = fovDegrees * (atan(1) * 4) / 180; //convert fov to radians
	for (int k = 0; k < sphereList->size(); k++) {

		Sphere currSphere = sphereList->at(k);

		currSphere.intersect();

		int i;
		int j;
		int h;
		int w;

		Vector origin(0, 0, 100);
		Vector d(0, 0, -1);
		Vector cs = currSphere.cs; //centre of the sphere
		double r = currSphere.r; //sphere radius
		Vector intersection;
		double t; //we find this
		Vector v(0, 0, 0);
		double a, b, c; //values for the quadratic equation
		Vector Lnorm; //light normal
		Vector light(400, -500, 1000); //light direction vector
		Vector snorm; //surface normal
		Vector reflection;

		const int phong = 64; //phong reflection coefficient

		double rCol;
		double gCol;
		double bCol;

		double colDiffuse;
		double colPhong;
		double colAmbient;

		double lightColRed;
		double lightColGreen;
		double lightColBlue;

		double ambColRed;
		double ambColGreen;
		double ambColBlue;

		double diffuseColRed;
		double diffuseColGreen;
		double diffuseColBlue;

		for (j = 0; j < hRes; j++) {
			for (i = 0; i < vRes; i++) {

				double x = 1 * (i - 0.5 * (imageWidth - 1));
				double y = 1 * (j - 0.5 * (imageHeight - 1));

				origin = Vector(x, y, 1);
				//origin.normalise();
				d.normalise();
				//origin = Vector(i - 500, j - 500, 1); //position of the pixel in 3D space, subtract 250 as 500 by 500 image
				v = origin.sub(cs);
				a = d.dot(d); //a = d*d
				b = 2 * (v.dot(d)); //b = 2*v*d
				c = v.dot(v) - r * r; //c = v*v - r*r
				double disc = b * b - 4 * a * c; //discriminant calculation
				if (disc < 0) {
					rCol = 0;
					gCol = 0;
					bCol = 0;
				}//if negative, we miss so no solution
				else {
					//update all values each loop
					lightColRed = 1;
					lightColGreen = 1;
					lightColBlue = 1;

					ambColRed = currSphere.colR;
					ambColGreen = currSphere.colG;
					ambColBlue = currSphere.colB;

					diffuseColRed = currSphere.colR;
					diffuseColGreen = currSphere.colG;
					diffuseColBlue = currSphere.colB;

					colPhong = .7;
					colDiffuse = .7;
					colAmbient = .3;

					ambColRed = ambColRed * colAmbient; //ambient is the colour intensity of each channel multiplied by the ambient coefficient
					ambColGreen = ambColGreen * colAmbient;
					ambColBlue = ambColBlue * colAmbient;

					//we have hit the circle
					t = (-b - sqrt(disc)) / 2 * a;
					intersection = origin.add(d.mul(t)); //calculate the point of intersection for this pixel

					//diffuse shading
					Lnorm = light.sub(intersection); //light normal calculated by light vector subtract intersection vector
					snorm = intersection.sub(cs); //surface normal is intersection vector subtract the centre of the circle

					Lnorm.normalise(); //normalise the light normal
					snorm.normalise(); //normalise the surface normal

					colDiffuse = snorm.dot(Lnorm) * colDiffuse; //colour intensity is the costheta, calculated by the dot product of surface normal and light normal

					if (colDiffuse <= 0) colDiffuse = 0; //if costheta is less that 0, it is in shadow
					if (colDiffuse > 1) colDiffuse = 1; //if not, use the costheta value

					diffuseColRed = diffuseColRed * colDiffuse; //diffuse is the amount reflected multiplied by the diffuse coefficient in each channel
					diffuseColGreen = diffuseColGreen * colDiffuse;
					diffuseColBlue = diffuseColBlue * colDiffuse;

					//add specular to colour here

					if (colDiffuse != 0) { //stop specular effect on the shadowed side of the sphere
						reflection = snorm.mul(2).mul(snorm.dot(Lnorm)).sub(Lnorm);//R = 2(N . L)N - L
					}
					else {
						colPhong = 0; //if in shadow, set specular component to 0
					}
					//add the colour to cosphy to the power of the phong coefficient
					intersection.normalise();
					reflection.normalise();

					colPhong = (colPhong * abs(pow(intersection.dot(reflection), phong)));

					rCol = ambColRed + diffuseColRed + colPhong * lightColRed; //total colour in each channel is the ambient + diffuse + phong, all multiplied by the light colour
					gCol = ambColGreen + diffuseColGreen + colPhong * lightColGreen;
					bCol = ambColBlue + diffuseColBlue + colPhong * lightColBlue;

					if (rCol > 1) rCol = 1;
					if (rCol < 0) rCol = 0;
					if (gCol > 1) gCol = 1;
					if (gCol < 0) gCol = 0;
					if (bCol > 1) bCol = 1;
					if (bCol < 0) bCol = 0;
					graphic.SetRGB(i, j, rCol * 255, gCol * 255, bCol * 255); //set the colour of the pixels depending on if sphere is intersected

				}

			}
		}
	}
}

cMain::cMain() : wxFrame(nullptr, wxID_ANY, "Testing Testing", wxPoint(hRes, vRes), wxSize(hRes, vRes))
{

	openConsole();

	graphic = new wxImage(hRes, vRes, true);

	//start rendering methods here...

	//startTriangleRayTrace(graphic, "dragon_vrip_res2.ply");
	startSphereRayTrace(graphic);
	startBVHRayTrace(graphic, "bun_zipper_res2.ply");
	
	Center();
}

//A method to start the BVH ray tracing
void cMain::startBVHRayTrace(wxImage* graphic, std::string fileName) {

	std::vector<Triangle> triangleList = plyRead(fileName).readFile();
	//create a base box for a root node of the hierarchy
	AABB baseBox = AABBintersect(triangleList);

	//create a node from the basebox
	BVHNode node = BVHNode(baseBox);
	node.setDepth(0);
	//create a tree from the root node
	BVHTree bvhtree = BVHTree(&node);
	std::vector<AABB> boxList;
	//add the base box to the box list
	boxList.push_back(baseBox);

	//print the number of triangles to render
	std::cout << "Base box triangle number: " << baseBox.getTriangleList().size() << std::endl;

	//construct the tree for the bounding volume hierarchy
	recursiveConstructBVH(&bvhtree, true);

	//print the depth of the tree
	std::cout << "Tree Depth: " << bvhtree.getRoot()->maxDepth(bvhtree.getRoot()) << std::endl;

	//print the depth of the right side of the tree
	std::cout << "Right Depth : " << bvhtree.getRoot()->getRight()->maxDepth(bvhtree.getRoot()->getRight()) << std::endl;

	//print the depth of the left side of the tree
	std::cout << "Left Depth : " << bvhtree.getRoot()->getLeft()->maxDepth(bvhtree.getRoot()->getLeft()) << std::endl;
	//add the leaves of the tree to a vector
	populateBoxList(&node);
	//print out the number of leaves in the tree
	std::cout << "Leaf number : " << leafBoxList.size();

	for (int i = 0; i < leafBoxList.size(); i++) {
		if (leafBoxList[i].getTriangleList().size() == 0) {
			std::cout << "what the heck" << '\n';
		}
		//std::cout << leafBoxList[i].getTriangleList().size() << '\n';
	}

	RayTraceBVH(graphic, triangleList, boxList, bvhtree);
}

//a method to start the sphere ray tracing
void cMain::startSphereRayTrace(wxImage* graphic) {

	Sphere sphere1 = Sphere(100, Vector(0, 0, 2), 1, 0, 0); //defines a red ball with radius 100 at the centre of the scene.
	Sphere sphere2 = Sphere(100, Vector(200, -100, -200), 0, 1, 0);
	Sphere sphere3 = Sphere(60, Vector(-200, -200, -300), 0, 0, 1);
	Sphere sphere4 = Sphere(80, Vector(-100, 150, -40), 1, 1, 0);
	Sphere sphere5 = Sphere(15, Vector(40, 20, 20), 0, 1, 1);
	Sphere sphere6 = Sphere(20, Vector(100, -30, -50), 1, 0, 1);
	Sphere sphere7 = Sphere(40, Vector(-210, -60, -40), .5, 1, 0);


	std::vector<Sphere> sphereList;
	sphereList.push_back(sphere5);
	sphereList.push_back(sphere1);
	sphereList.push_back(sphere4);
	sphereList.push_back(sphere7);
	sphereList.push_back(sphere6);
	sphereList.push_back(sphere2);
	sphereList.push_back(sphere3);


	RayTraceSphere(&sphereList, *graphic);
}

//a method to start triangle ray tracing (no acceleration structure)
void cMain::startTriangleRayTrace(wxImage* graphic, std::string fileName) {
	std::vector<Triangle> triangleList = plyRead(fileName).readFile();

	std::cout << "Triangle mesh size: " << triangleList.size() << std::endl;

	RayTraceTriangleMesh(graphic, triangleList);
}

cMain::~cMain() 
{
}

void cMain::OnPaint(wxPaintEvent &evt)
{
	wxPaintDC dc(this);
	bmp = new wxBitmap(*graphic);
	dc.DrawBitmap(*bmp, 0, 0, false);
	//fix memory leaks here
	delete bmp;
	delete graphic;
	evt.Skip();
}


void cMain::OnButtonClicked(wxCommandEvent &evt)
{
	//m_list1->AppendString(m_txt1->GetValue());
	
	evt.Skip();
}

//Method to open the console using a non-windowed application
void cMain::openConsole() {
	AllocConsole();
	FILE* pNewStdout = nullptr;
	FILE* pNewStderr = nullptr;
	FILE* pNewStdin = nullptr;

	::freopen_s(&pNewStdout, "CONOUT$", "w", stdout);
	::freopen_s(&pNewStderr, "CONOUT$", "w", stderr);
	::freopen_s(&pNewStdin, "CONIN$", "r", stdin);

	// Clear the error state for all of the C++ standard streams. Attempting to accessing the streams before they refer
	// to a valid target causes the stream to enter an error state. Clearing the error state will fix this problem,
	// which seems to occur in newer version of Visual Studio even when the console has not been read from or written
	// to yet.

	std::cout.clear();
	std::cerr.clear();
	std::cin.clear();

	std::wcout.clear();
	std::wcerr.clear();

	std::wcin.clear();
}

