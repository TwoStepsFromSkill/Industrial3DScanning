#include "KDTree.h"
#include "Matrix.h"
#include "SVD.h"

#include <algorithm>
#include <limits>
#include <cmath>
#include <iostream>

//Constructor
Node::Node()
	: median(0)
	, leftChild(nullptr)
	, rightChild(nullptr)
	, ptrFirstPoint(nullptr)
	, ptrLastPoint(nullptr)
{}

//Deconstructor
Node::~Node()
{
	if (leftChild)
		delete leftChild;
	if (rightChild)
		delete rightChild;
}

//Build Functions
Node* KDTree::buildKDTree(Point3d* rangeBegin, Point3d* rangeEnd, unsigned int depth)
{
	//compute how many points there are in the given range
	unsigned int numPoints = (rangeEnd - rangeBegin); //Pointers can be subtracted
													  //after sorting the median position is at the half of the range
	unsigned int medianPosition = numPoints / 2;
	//the median value in the current range
	double medianValue = 0; //initialize with 0
							//in depth=0 we are processing x-values, at depth=1 y-values, at depth=2 z-Values,
							//at depth=3 x-values ... and so on
	unsigned int currentDimension = depth % 3; //Modulo-Operation
	Point3d* pointAtMedianPosition;
	//find median depending on the current tree depth
	if (currentDimension == 0)
	{
		std::sort(rangeBegin, rangeEnd, sortByXvalue);
		pointAtMedianPosition = (rangeBegin + medianPosition);
		medianValue = pointAtMedianPosition->x;
	}
	else if (currentDimension == 1)
	{
		std::sort(rangeBegin, rangeEnd, sortByYvalue);
		pointAtMedianPosition = (rangeBegin + medianPosition);
		medianValue = pointAtMedianPosition->y;
	}
	else
	{
		std::sort(rangeBegin, rangeEnd, sortByZvalue);
		pointAtMedianPosition = (rangeBegin + medianPosition);
		medianValue = pointAtMedianPosition->z;
	}

	//OK range is sorted, we create a new childNode
	Node* childNode = new Node; //create new node
	childNode->median = medianValue;
	childNode->ptrFirstPoint = rangeBegin;
	childNode->ptrLastPoint = rangeEnd;

	if (numPoints > 1)
	{
		if (numPoints % 2 != 0)
		{
			pointAtMedianPosition++;
		}

		childNode->leftChild = buildKDTree(rangeBegin, pointAtMedianPosition, depth + 1);
		childNode->rightChild = buildKDTree(pointAtMedianPosition, rangeEnd, depth + 1);
	}

	return childNode;
}

bool KDTree::sortByXvalue(const Point3d& p1, const Point3d& p2)
{
	if (p1.x < p2.x)
		return true;
	else
		return false;
}
bool KDTree::sortByYvalue(const Point3d& p1, const Point3d& p2)
{
	if (p1.y < p2.y)
		return true;
	else
		return false;
}
bool KDTree::sortByZvalue(const Point3d& p1, const Point3d& p2)
{
	if (p1.z < p2.z)
		return true;
	else
		return false;
}

//Range Query
std::vector<Point3d> queryRange(Node* tree, const double minMax[6])
{
	std::vector<Point3d> result;

	queryRange_impl(tree, minMax, 0, result);

	return result;
}

/**
@brief	initializes the spherical range query
@param	tree is the K-D-Tree which stores the point cloud
@param	center is the middle of the sphere
@param	radius specifies the search area around the center point
@param	flagging is false by default and only used for efficient thinning
*/
std::vector<Point3d*> queryRangeSphere(Node* tree, Point3d* center, const double radius, bool flagging)
{
	std::vector<Point3d*> result;

	queryRangeSphere_impl(tree, center, radius, 0, result, flagging);

	return result;
}

void queryRange_impl(Node* tree, const double minMax[6], unsigned int depth, std::vector<Point3d>& out)
{
	// Stepped to deep
	if (!tree)
		return;

	// Reached a leaf -> check if points are inside range and add the ones that are
	if (!tree->leftChild && !tree->rightChild)
	{
		auto* start = tree->ptrFirstPoint;

		while (start != tree->ptrLastPoint)
		{
			Point3d pt = *start;

			if (pt.x > minMax[0] && pt.x < minMax[1]
				&& pt.y > minMax[2] && pt.y < minMax[3]
				&& pt.z > minMax[4] && pt.z < minMax[5])
			{
				out.push_back(*start);
			}

			++start;
		}
	}
	else
	{
		unsigned int dimension = depth % 3;

		if (minMax[dimension * 2] <= tree->median)
			queryRange_impl(tree->leftChild, minMax, dimension + 1, out);
		if (minMax[dimension * 2 + 1] > tree->median)
			queryRange_impl(tree->rightChild, minMax, dimension + 1, out);
	}
}

/**
@brief	implements the spherical range query
@param	tree is the K-D-Tree which stores the point cloud
@param	center is the middle of the sphere
@param	radius specifies the search area around the center point
@param	depth stores the current level of the traversed tree
@param	out stores the points founds in the neighborhood
@param	flagging is false by default and only used for efficient thinning
*/
void queryRangeSphere_impl(Node* tree, Point3d* center, const double radius, unsigned int depth, std::vector<Point3d*>& out, bool flagging)
{
	// Stepped to deep
	if (!tree)
		return;

	// Reached a leaf -> check if points are inside the spherical neighborhood and add the ones that are
	if (!tree->leftChild && !tree->rightChild)
	{
		auto* start = tree->ptrFirstPoint;
		
		//iterate once
		while (start != tree->ptrLastPoint)
		{
			Point3d pt = *start;
			
			//calculate distance between pt and center
			const double dist = distance3d(*center, pt);

			//flag mode on (required for efficient thinning)
			if (flagging)
			{
				if (dist <= radius)
				{
					//add leaf to output vector if it lies within sphere(center, radius)
					out.push_back(tree->ptrFirstPoint);

					//flag detected point if it's not the current center
					if (dist != 0)
					{
						tree->ptrFirstPoint->flag_ignore = true;
					}
				}
			}
			else //flag mode off
			{
				//add leaf to output vector if it lies within sphere(center, radius)
				if (dist <= radius)
				{
					out.push_back(tree->ptrFirstPoint);
				}
			}

			++start;
		}
	}
	else //traverse tree
	{
		unsigned int dimension = depth % 3;

		//switch dimension for comparison with median
		double cmpr;

		switch (dimension)
		{
			case 0: cmpr = center->x; break;
			case 1: cmpr = center->y; break;
			case 2: cmpr = center->z; break;

			default: break;
		}

		if (cmpr - radius <= tree->median)
			queryRangeSphere_impl(tree->leftChild, center, radius, dimension + 1, out, flagging);
		if (cmpr + radius > tree->median)
			queryRangeSphere_impl(tree->rightChild, center, radius, dimension + 1, out, flagging);
	}
}

//Nearest Neighbor Search
Point3d nearestNeighbor_daniel(Node* tree, const Point3d& queryPoint)
{
	//create "infinity" to use as first reference
	double minDist = std::numeric_limits<double>::max();
	
	Point3d minPoint(std::numeric_limits<double>::lowest(),
		std::numeric_limits<double>::lowest(),
		std::numeric_limits<double>::lowest());
	
	//create array from query point data for easy access
	double query[3] = { queryPoint.x, queryPoint.y, queryPoint.z };

	//call search function
	nearestNeighbor_daniel_impl(tree, query, &minDist, &minPoint, 0);

	return minPoint;
}

void nearestNeighbor_daniel_impl(Node* tree, const double* queryPoint, double* minDist, Point3d* minPoint, unsigned int depth)
{
	if (!tree)
		return;

	if (!tree->leftChild && !tree->rightChild)
	{
		auto* start = tree->ptrFirstPoint;

		while (start != tree->ptrLastPoint)
		{
			Point3d pt = *start;
			double dist = std::sqrt(sqr(pt.x - queryPoint[0])
				+ sqr(pt.y - queryPoint[1])
				+ sqr(pt.z - queryPoint[2]));
			if (dist < *minDist)
			{
				*minDist = dist;
				*minPoint = pt;
			}

			++start;
		}
	}
	else
	{
		unsigned int dim = depth % 3;

		if (queryPoint[dim] <= tree->median)
		{
			nearestNeighbor_daniel_impl(tree->leftChild, queryPoint, minDist, minPoint, depth + 1);

			if ((queryPoint[dim] + *minDist) > tree->median)
				nearestNeighbor_daniel_impl(tree->rightChild, queryPoint, minDist, minPoint, depth + 1);
		}
		else
		{
			nearestNeighbor_daniel_impl(tree->rightChild, queryPoint, minDist, minPoint, depth + 1);

			if ((queryPoint[dim] - *minDist) <= tree->median)
				nearestNeighbor_daniel_impl(tree->leftChild, queryPoint, minDist, minPoint, depth + 1);
		}
	}
}

//Smoothing
std::vector<Point3d> computeSmoothing(std::vector<Point3d> points, Node* tree, const double radius)
{
	/*double radius;
	m_smoothingWidget->getRadius(&radius);*/

	std::vector<Point3d> smoothedPoints;

	bool useGaussSmoothing = false;

	if (useGaussSmoothing)
	{
		smoothedPoints = smoothPointsGaussian(points, tree, radius);
	}
	else
	{
		smoothedPoints = smoothPointsAverage(points, tree, radius);
	}

	//std::vector<unsigned char> colors = computeColorsForSmoothing(smoothedPoints);

	/*m_glWidget->setSmoothedPoints(smoothedPoints);
	m_glWidget->setPointColors(colors);*/

	/*emit drawingTemporaryChanged(false);
	emit drawingMainPointCloudChanged(false);
	emit drawingSmoothedPointsChange(true);*/

	return smoothedPoints;
}

std::vector<Point3d> smoothPointsAverage(std::vector<Point3d>& points, Node* rootNode, double radius)
{
	std::vector<Point3d> smoothedPoints;
	smoothedPoints.resize(points.size());

	#pragma omp parallel for
	for (int i = 0; i < points.size(); ++i)
	{
		Point3d& point = points[i];

		Point3d smoothedPointAv;
		std::vector<Point3d*> neighborPoints = queryRangeSphere(rootNode, &point, radius, false);

		for (std::size_t j = 0; j < neighborPoints.size(); ++j)
		{
			smoothedPointAv += *neighborPoints[j];
		}

		smoothedPointAv *= 1.0 / neighborPoints.size();
		smoothedPoints[i] = smoothedPointAv;
	}

	return smoothedPoints;
}

/**
@brief	smoothPointsGaussian is a function that smoothes the points with an gaussian kernel
		in a spherical neighborhood
@param	points is a vector which contains all given points
@param	rootNode is a 3dTree which contains all given points in a sorted order
@param	radius defines the range for the spherical neighboorhood search
@return is a vector which contains all smoothed points
*/
std::vector<Point3d> smoothPointsGaussian(std::vector<Point3d>& points, Node* rootNode, double radius)
{
	std::vector<Point3d> smoothedPoints;
	smoothedPoints.resize(points.size());

	#pragma omp parallel for
	for (int i = 0; i < points.size(); ++i)
	{
		Point3d& point = points[i];

		Point3d smoothedPointSum;
		std::vector<Point3d*> neighborPoints = queryRangeSphere(rootNode, &point, radius, false);
		double sumWeights = 0;

		for (std::size_t j = 0; j < neighborPoints.size(); ++j)
		{
			double distance = sqDistance3d(point, *neighborPoints[j]);
			double weight = std::exp((-distance) / radius);
			smoothedPointSum += (*neighborPoints[j] * weight);
			sumWeights += weight;
		}
		if (sumWeights != 0)
			smoothedPoints[i] = smoothedPointSum * (1 / sumWeights);
	}

	return smoothedPoints;
}

/*
std::vector<unsigned char> computeColorsForSmoothing(std::vector<Point3d> points, const std::vector<Point3d>& smoothedPoints)
{
	std::vector<unsigned char> colors(3 * points.size(), 0);
	std::vector<double> distances(points.size());

	#pragma omp parallel for
	for (int i = 0; i < points.size(); ++i)
	{
		distances[i] = sqDistance3d(points[i], smoothedPoints[i]);
	}

	auto minmax = std::minmax_element(distances.begin(), distances.end());
	double min = *minmax.first;
	double max = *minmax.second;

	// Grey scale
	// #pragma omp parallel for
	//     for (int i = 0; i < distances.size(); ++i)
	//     {
	//         unsigned char factor = ((distances[i] - min) / (max - min)) * 255;
	//         colors[i*3] = factor;
	//         colors[i*3 + 1] = 0;
	//         colors[i*3 + 2] = 0;
	//     }

	// Between two colors
	//     constexpr unsigned char from[3] {31, 28, 24};
	//     constexpr unsigned char to[3] {142, 14, 0};
	//
	// #pragma omp parallel for
	//     for (int i = 0; i < distances.size(); ++i)
	//     {
	//         double factor = ((distances[i] - min) / (max - min));
	//         colors[i*3] = (1.0 - factor)*from[0] + factor*to[0];
	//         colors[i*3 + 1] = (1.0 - factor)*from[1] + factor*to[1];
	//         colors[i*3 + 2] = (1.0 - factor)*from[2] + factor*to[2];
	//     }

	// Multiple colors (equally spaced)

	// Rainbow
	//     const std::vector<std::array<unsigned char, 3>> lookup {{0,0,255},
	//                                                             {0,255,0},
	//                                                             {255,255,0},
	//                                                             {255,127,0},
	//                                                             {255,0,0}
	//     };
	// Heat
	const std::vector<std::array<unsigned char, 3>> lookup{ { 0,0,0 },
	{ 122,0,0 },
	{ 255,229,33 },
	{ 255,255,255 }
	};
	const double spacing = 1.0 / (lookup.size() - 1);

	#pragma omp parallel for
	for (int i = 0; i < distances.size(); ++i)
	{
		double factor = ((distances[i] - min) / (max - min));

		std::size_t startIdx = std::floor((lookup.size() - 1) * factor);
		std::size_t endIdx = startIdx + 1;

		double startFactor = startIdx * spacing;
		double endFactor = endIdx * spacing;

		double factorBetween = (factor - startFactor) / (endFactor - startFactor);

		colors[i * 3] = (1.0 - factorBetween)*lookup[startIdx][0] + factorBetween*lookup[endIdx][0];
		colors[i * 3 + 1] = (1.0 - factorBetween)*lookup[startIdx][1] + factorBetween*lookup[endIdx][1];
		colors[i * 3 + 2] = (1.0 - factorBetween)*lookup[startIdx][2] + factorBetween*lookup[endIdx][2];
	}

	return colors;
}*/

std::vector<double> colorFromGradientHSV(double index)
{
	if (index < 0) index = 0;
	if (index > 1) index = 1;

	const double H(240 * (1 - index));
	const int hi(std::floor(H / 60));
	const double f(H / 60 - hi);
	const double V(1);
	const double S(1);
	const double p(V*(1 - S));
	const double q(V*(1 - S*f));
	const double t(V*(1 - S*(1 - f)));

	std::vector<double> color;

	switch (hi)
	{
	case 1:
		color = { 255 * q, 255 * V, 255 * p }; break;
	case 2:
		color = { 255 * p, 255 * V, 255 * t }; break;
	case 3:
		color = { 255 * p, 255 * q, 255 * V }; break;
	case 4:
		color = { 255 * t, 255 * p, 255 * V }; break;
	case 5:
		color = { 255 * V, 255 * p, 255 * q }; break;

	default:
		color = { 255 * V, 255 * t, 255 * p }; break;
	}

	return color;
}

//Thinning
/**
@brief	Thins a given pointcloud (stored in K-D-Tree) by a given radius.
@param	globalTree root node of the K-D-Tree that stores the original point data
@param	subTree node of the currently (used for recursive traversion)
@param	radius specifies the degree of thinning
@param	output stores the remaining points
*/
void homogeneousThinning(Node* globalTree, Node* subTree, const double radius, std::vector<Point3d>& output) //needs to distinguish between global and sub tree
{
	//stepped to deep
	if (!subTree)
		return;

	if (!subTree->leftChild && !subTree->rightChild) //reached a leaf
	{
		if (!subTree->ptrFirstPoint->flag_ignore) //leaf hasn't been flagged
		{
			//store leaf value as point
			Point3d curPoint = *subTree->ptrFirstPoint;

			//add point to output vector
			output.push_back(curPoint);

			//get points in spherical neighborhood(curPoint, radius) and flag them
			std::vector<Point3d*> neighborhood = queryRangeSphere(globalTree, &curPoint, radius, true);
			//std::cout << "\nneighborhood.size: " << neighborhood.size();
		}
	}
	else //traverse tree
	{
		homogeneousThinning(globalTree, subTree->leftChild, radius, output);
		homogeneousThinning(globalTree, subTree->rightChild, radius, output);
	}
}

//Best-Fit Plane
/**
@brief	Computes centroid p of all given points.
@param	given point cloud
@return	a Point3d = centroid of point cloud
*/
Point3d getCentroid(std::vector<Point3d> points)
{
	Point3d result;

	const int N = points.size();

	for each (Point3d pt in points)
	{
		result.x += pt.x;
		result.y += pt.y;
		result.z += pt.z;
	}

	result *= 1 / N;

	return result;
}

/**
@brief	Calculates the covariance matrix of a given point cloud with a specified centroid.
@param	given point cloud
@param	centroid of given point cloud
@return	covariance matrix of given point cloud
*/
Matrix getCovarianceMatrix(std::vector<Point3d> points, Point3d centroid)
{
	const int N = points.size();

	double x_, y_, z_;

	double	Mxx = 0, Mxy = 0, Mzx = 0, 
			Myy = 0, Myz = 0,
			Mzz = 0;

	for (size_t i = 0; i < N; i++)
	{
		x_ = points[i].x - centroid.x;
		y_ = points[i].y - centroid.y;
		z_ = points[i].z - centroid.z;

		Mxx += x_ * x_;
		Mxy += x_ * y_;
		Mzx += z_ * x_;
		Myy += y_ * y_;
		Myz += y_ * z_;
		Mzz += z_ * z_;
	}

	//create 3x3 output matrix
	Matrix result = Matrix(3, 3);

	//fill matrix
	result(0, 0) = Mxx; result(0, 1) = Mxy; result(0, 2) = Mzx;
	result(1, 0) = Mxy; result(1, 1) = Myy; result(1, 2) = Myz;
	result(2, 0) = Mzx; result(2, 1) = Myz; result(2, 2) = Mzz;

	return result;
}

/**
@brief	Calculates the direction of the best-fit plane by SVD decomposition.
@param	covariance matrix of given point cloud
@return	direction of best-fit plane normal
*/
Point3d getDirectionOfBestFitPlane(Matrix covariance)
{
	std::vector<double> S;

	Matrix V;

	SVD::decomposeMatrix(covariance, S, V);

	Point3d result = Point3d(	covariance(0, 2),
								covariance(1, 2),
								covariance(2, 2));

	return result;
}
