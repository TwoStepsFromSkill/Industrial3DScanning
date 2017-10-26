#pragma once
#include <algorithm>

#include "Point3d.h"

struct Node
{
	
	double median;
	Node* leftChild;
	Node* rightChild;
	Point3d* ptrFirstPoint;
	Point3d* ptrLastPoint;
	~Node()
	{
		delete leftChild;
		delete rightChild;
	}
};

class KDTree
{
public:

	static Node* buildKDTree(Point3d* rangeBegin, Point3d* rangeEnd, unsigned intdepth);
	
	static bool sortByXvalue(const Point3d& p1, const Point3d& p2);
	
	static bool sortByYvalue(const Point3d& p1, const Point3d& p2);
	
	static bool sortByZvalue(const Point3d& p1, const Point3d& p2);
};

