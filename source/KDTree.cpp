#include "KDTree.h"

#include <algorithm>
#include <limits>
#include <cmath>

Node::Node()
    : median(0)
    , leftChild(nullptr)
    , rightChild(nullptr)
    , ptrFirstPoint(nullptr)
    , ptrLastPoint(nullptr)
{}

Node::~Node()
{
    if (leftChild)
        delete leftChild;
    if (rightChild)
        delete rightChild;
}

Point3d* givenPoint = new Point3d(0.0010,0.0347,-0.0236);

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
	else if(currentDimension == 1)
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

std::vector<Point3d> queryRange(Node* tree, const double minMax[6])
{
    std::vector<Point3d> result;
    queryRange_impl(tree, minMax, 0, result);

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

        if (minMax[dimension*2] <= tree->median)
            queryRange_impl(tree->leftChild, minMax, dimension + 1, out);
        if (minMax[dimension*2 + 1] > tree->median)
            queryRange_impl(tree->rightChild, minMax, dimension + 1, out);
    }
}

Point3d nearestNeighbor_daniel(Node* tree, const Point3d& queryPoint)
{
    double minDist = std::numeric_limits<double>::max();
    Point3d minPoint(std::numeric_limits<double>::lowest(),
                     std::numeric_limits<double>::lowest(),
                     std::numeric_limits<double>::lowest());
    double query[3] = {queryPoint.x, queryPoint.y, queryPoint.z};
	nearestNeighbor_daniel_impl(tree, query, &minDist, &minPoint, 0);
	return minPoint;
}
void nearestNeighbor_daniel_impl(Node* tree, const double* queryPoint, double* minDist,
	Point3d* minPoint, unsigned int depth)
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
Point3d findNearestPoint_Elke(Node* tree)
{	
	Point3d currentPoint = *tree->ptrFirstPoint;
	double currDist = distance3d(currentPoint, *givenPoint);
	//givenPoint = new Point3d(0.9, 0.9, 0.9);
	
	insertGivenPoint(tree, &currDist, 0, &currentPoint);
	findNearestPointRecursiveley_Elke(tree, &currDist,0,&currentPoint);

	return currentPoint;
}

void findNearestPointRecursiveley_Elke(Node* tree, double* currDist, unsigned int depth, Point3d* currentPoint)
{
	// Stepped to deep
	if (!tree)
		return;

	unsigned int dimension = depth % 3;
	double givenPointDimension = getDimension(*givenPoint, depth);
	double distToMedian = givenPointDimension - tree->median;

	// If leaf reached -> Check left child
	if (!tree->leftChild && !tree->rightChild)
	{		
		if (distance3d(*tree->ptrFirstPoint, *givenPoint) < *currDist)
		{
			*currDist = distance3d(*tree->ptrFirstPoint, *givenPoint);
			*currentPoint = *tree->ptrFirstPoint;
		}
		return;
	}
	
	// givenPoint is actually on the left side of the tree
	if (distToMedian <= 0)
	{
		findNearestPointRecursiveley_Elke(tree->leftChild, currDist, depth + 1, currentPoint);
		if (*currDist > abs(distToMedian))
		{
			findNearestPointRecursiveley_Elke(tree->rightChild, currDist, depth + 1, currentPoint);
		}
	}
	else 
	{
		findNearestPointRecursiveley_Elke(tree->rightChild, currDist, depth + 1, currentPoint);
		if (*currDist > abs(distToMedian))
		{
			findNearestPointRecursiveley_Elke(tree->leftChild, currDist, depth + 1, currentPoint);
		}		
	}	
}

void insertGivenPoint(Node* tree, double* currDist, unsigned int depth, Point3d* currentPoint)
{
	if (!tree->leftChild && !tree->rightChild)
	{
		currentPoint = tree->ptrFirstPoint;
		*currDist = distance3d(*currentPoint, *givenPoint);
		return;
	}
	else
	{
		unsigned int dimension = depth % 3;
		double givenPointDimension = getDimension(*givenPoint, depth);
		if (tree->median <= givenPointDimension)
		{
			depth++;
			insertGivenPoint(tree->leftChild, currDist, depth, currentPoint);
		}
		else
		{
			depth++;
			insertGivenPoint(tree->rightChild, currDist, depth, currentPoint);
		}
	}
}

double getDimension(Point3d point, unsigned int depth)
{
	unsigned int dimension = depth % 3;
	double currentPosition = 0;
	switch (dimension)
	{
	case 0:
		currentPosition = point.x;
		break;
	case 1:
		currentPosition = point.y;
		break;
	case 2:
		currentPosition = point.z;
		break;
	default:
		break;
	}
	return currentPosition;
}

    

