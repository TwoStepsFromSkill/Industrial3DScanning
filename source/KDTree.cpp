#include "KDTree.h"

#include <algorithm>

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

    auto recursiveDown = [&minMax, &tree, &out](unsigned int dimension)
    {
        if (minMax[dimension*2] <= tree->median)
            queryRange_impl(tree->leftChild, minMax, dimension + 1, out);
        if (minMax[dimension*2 + 1] > tree->median)
            queryRange_impl(tree->rightChild, minMax, dimension + 1, out);
    };

    recursiveDown(depth % 3);
}