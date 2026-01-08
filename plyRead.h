//Sam Harry, Student Number 1901522

#pragma once
#include <iostream>
#include <fstream>
#include "Triangle.h"
#include <vector>

class plyRead {

public: 
	std::string filename;

public:

	plyRead(std::string filename);
	
	~plyRead();

	std::vector<Triangle> readFile();
};