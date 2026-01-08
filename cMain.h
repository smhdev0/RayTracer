//Sam Harry, Student Number 1901522

#pragma once
#include "wx/wx.h"
#include <vector>
#include "Vector.h"
#include "Sphere.h"
#include "Triangle.h"
#include "AABB.h"
#include "BVHNode.h"
#include "BVHTree.h"

class cMain : public wxFrame
{
public :
	cMain();
	~cMain();
	static void RayTraceSphere(std::vector<Sphere> *sphereList, wxImage &graphic);
	static void RayTraceTriangleMesh(wxImage* graphic, std::vector<Triangle>& triangleList);
	static void RayTraceBVH(wxImage* graphic, std::vector<Triangle>& triangleList, std::vector<AABB>& AABBlist, BVHTree& tree);
	void startBVHRayTrace(wxImage* graphic, std::string fileName);
	void startSphereRayTrace(wxImage* graphic);
	void startTriangleRayTrace(wxImage* graphic, std::string fileName);
	static bool isWithinAABB(AABB * aabb, double x, double y);
	void populateBoxList(BVHNode* node);
	static BVHNode * getNodeHit(BVHNode * node, double x, double y);
	std::vector<AABB> testSplit(std::vector<Triangle>& triangleList);
	std::vector<Triangle> getAABBTriangles(AABB * box, std::vector<Triangle> & triangleList);
	double getAverageTriangleCentroid(AABB* box, bool isVerticalSplit);
	void recursiveConstructBVHBox(BVHNode * node, BVHTree * tree, bool isVertical);
	void recursiveConstructBVH(BVHTree * tree, bool isVertical);
	static AABB AABBintersect(std::vector<Triangle> & triangleList);
	static std::vector<AABB> constructBVH(std::vector<Triangle> &triangleList, wxImage & graphic);
	void openConsole();

public :
	wxButton *m_btn1 = nullptr;
	wxTextCtrl *m_txt1 = nullptr;
	wxListBox *m_list1 = nullptr;
	wxImage *graphic = nullptr;
	wxBitmap *bmp = nullptr;
	wxWindow *window = nullptr;
	wxPaintDC *dc = nullptr;

	void OnPaint(wxPaintEvent & evt);
	void OnButtonClicked(wxCommandEvent &evt);

	wxDECLARE_EVENT_TABLE();
};

