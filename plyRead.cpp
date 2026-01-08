//Sam Harry, Student Number 1901522

#include <iostream>
#include <cmath>
#include <fstream>
#include "plyread.h"
#include <string>
#include <sstream>
#include "Triangle.h"
#include <vector>
#include "Point.h"
#include <istream>



plyRead::plyRead(std::string filenamein)
	: filename(filenamein) {}

plyRead::~plyRead() {

}

std::vector<Triangle> plyRead::readFile() {

	//
	double x, y, z;
	int a, b, c;
	std::string line;
	std::ifstream file(this->filename);
	int count = 0;
	std::vector<Point> vertexList;
	std::vector<Triangle> triangleList;
	Triangle triangle;

	while (getline(file, line)) {
		if (line[0] == '3' && line[1] ==' ') {
			
				std::istringstream iss(line);
				iss.ignore();
				iss >> a;
				iss.ignore();
				iss >> b;
				iss.ignore();
				iss >> c;
				
				triangle = Triangle(vertexList[a], vertexList[b], vertexList[c]);
				triangleList.push_back(triangle);
			
		}
		else {
			std::istringstream iss(line);
			iss >> x;
			iss.ignore();
			iss >> y;
			iss.ignore();
			iss >> z;
			Point a = Point(x, y, z);

			vertexList.push_back(a);

		}
		count++;

		//std::cout << line << '\n';
	}

	



	file.close();
	return triangleList;

}