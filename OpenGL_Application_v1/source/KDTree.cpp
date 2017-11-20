#include "KDTree.h"

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
			//std::cout << "\ncenter: " << Point3dToString(center);
			//std::cout << "\n    pt: " << Point3dToString(&pt);
			
			//calculate distance between pt and center
			const double dist = distance3d(*center, pt);
			//std::cout << "\ndist = " << dist;

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
						//std::cout << "\nflagged: " << Point3dToString(tree->ptrFirstPoint);
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

//Thinning
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
