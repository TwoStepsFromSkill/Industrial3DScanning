#ifndef KD_TREE_H
#define KD_TREE_H

#include "Point3d.h"

#include <vector>

struct Node
{
	Node();
	~Node();

	double median;
	Node* leftChild;
	Node* rightChild;
	Point3d* ptrFirstPoint;
	Point3d* ptrLastPoint;
};

struct KDTree
{
	static Node* buildKDTree(Point3d* rangeBegin, Point3d* rangeEnd, unsigned int depth);

	static bool sortByXvalue(const Point3d& p1, const Point3d& p2);
	static bool sortByYvalue(const Point3d& p1, const Point3d& p2);
	static bool sortByZvalue(const Point3d& p1, const Point3d& p2);
};

//Boxed Range Query
std::vector<Point3d> queryRange(Node* tree, const double minMax[6]);
void queryRange_impl(Node* tree, const double minMax[6], unsigned int depth, std::vector<Point3d>& out);

//Spherical Range Query
std::vector<Point3d*> queryRangeSphere(Node* tree, Point3d* center, const double radius, bool flagging);
void queryRangeSphere_impl(Node* tree, Point3d* center, const double radius, unsigned int depth, std::vector<Point3d*>& out, bool flagging);

//Nearest Neighbor Search
Point3d nearestNeighbor_daniel(Node* tree, const Point3d& queryPoint);
void nearestNeighbor_daniel_impl(Node* tree, const double* queryPoint, double* minDist, Point3d* minPoint, unsigned int depth);

//Thinning
void homogeneousThinning(Node* globalTree, Node* subTree, const double radius, std::vector<Point3d>& output);

#endif // KD_TREE_H