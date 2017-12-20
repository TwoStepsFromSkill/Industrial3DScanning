#ifndef KD_TREE_H
#define KD_TREE_H

#include "Point3d.h"
#include "Matrix.h"

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

//Smoothing
std::vector<Point3d> computeSmoothing(std::vector<Point3d> points, Node* tree, const double radius);

std::vector<Point3d> smoothPointsAverage(std::vector<Point3d>& points, Node* rootNode, double radius);
std::vector<Point3d> smoothPointsGaussian(std::vector<Point3d>& points, Node* rootNode, double radius);

//std::vector<unsigned char> computeColorsForSmoothing(const std::vector<Point3d>& smoothedPoints);
std::vector<double> colorFromGradientHSV(double index);

//Thinning
/**
@image html Thinning_comparison.jpg
*/
void homogeneousThinning(Node* globalTree, Node* subTree, const double radius, std::vector<Point3d>& output);

//Best-Fit Plane
Point3d getCentroid(std::vector<Point3d>* points);

Matrix getCovarianceMatrix(std::vector<Point3d>* points, Point3d centroid);

std::vector<Point3d> getDirectionsOfBestFit(Matrix covariance, std::vector<Point3d>* eigenVectorsOutput);

std::vector<Point3d> getCornersOfBestFit(std::vector<Point3d>* points, Point3d centroid, std::vector<Point3d>* eigenVectors);

//Best-Fit Sphere
double getMeanDisToPoint(std::vector<Point3d>* points, Point3d point);

std::vector<double> getDistancesForJacobi(std::vector<Point3d>* points, Point3d center, double radius);

Matrix getJacobiMatrix(std::vector<Point3d>* points, Point3d center, std::vector<double> distances);

double getVariance(const std::vector<double>* data);

#endif // KD_TREE_H