#ifndef KD_TREE_H
#define KD_TREE_H

#include "Point3d.h"

#include <vector>

/**
 * @brief Node of a KDTree.
 */
struct Node
{
    /**
     * @brief Default constructor.
     */
	Node();
    /**
     * @brief Destructor.
     */
    ~Node();

	double median;
	Node* leftChild;
	Node* rightChild;
	Point3d* ptrFirstPoint;
	Point3d* ptrLastPoint;
};

struct KDTree
{
    /**
	@brief function which creats an 3DTree
	@details on the left side of the tree are the points with small x,y,z values
	@param rangeBegin is a pointer which is pointing on an Point3d
	@param rangeEnd is a pointer which is pointing on the adress behind the last Point3d
	@param depth is the depth of the tree with which the actual dimension can be find out
	@return the root node of the 3DTree
    */
	static Node* buildKDTree(Point3d* rangeBegin, Point3d* rangeEnd, unsigned int depth);

    /**
	@brief sorts two Point3d by their x value
    */
	static bool sortByXvalue(const Point3d& p1, const Point3d& p2);
    /**
    @brief sorts two Point3d by their y value
    */
	static bool sortByYvalue(const Point3d& p1, const Point3d& p2);
    /**
    @brief sorts two Point3d by their z value
    */
	static bool sortByZvalue(const Point3d& p1, const Point3d& p2);
};

/** @defgroup group_kdtree KDTree methods
 *  @{
 */

/**
 * @brief Computes all points that lie in the specified range.
 * @param[in] tree Root node of the KDTree.
 * @param[in] minMax 6-element array containing the min and max coordiantes of the range in x, y and
 *                      and z axis. Storage layout [xMin, xMax, yMin, yMax, zMin, zMax]
 * @return Returns a new vector with copies of the point positions.
 */
std::vector<Point3d> queryRange(Node* tree, const double minMax[6]);
/**
 * @brief Actual implementation of queryRange function.
 * @param[in] tree Current node in the KDTree.
 * @param[in] minMax 6-element array containing the min and max coordiantes of the range.
 * @param[in] depth Depth of the current node in the tree (needed for coordiante selection)
 * @param[in,out] out Vector containing the positions of points inside the range. Might be adapted
 *                  by this function if tree is a leaf and points are inside the range.
 */
void queryRange_impl(Node* tree, const double minMax[6], unsigned int depth,
                     std::vector<Point3d>& out);

/**
 * @brief function which calls the recursiv function that calculate all neighboorpoints of the centerPoint
 * @param tree is a pointer on the root node of a 3DTree
 * @param radius is the radius in which the function search spherical for the neighboorpoints of the centerPoint
 * @param centerPoint is the Point3d for which the neighboorpoints should be found
 * @return a vector of Point3d which contains all neighboorpoints
 */
std::vector<Point3d> queryRadius(Node* tree, const double radius, Point3d centerPoint);
void queryRadius_impl(Node* tree, const double radius, unsigned int depth,
	std::vector<Point3d>& out, Point3d centerPoint);

/**
 * @brief Find nearest neighbor of point in a KDTree.
 * @param[in] tree Root node of the KDTree.
 * @param[in] queryPoint Find nearest neighbor of this point.
 * @return Returns a copy of the position of the nearest neighbor of queryPoint.
 */
Point3d nearestNeighbor_daniel(Node* tree, const Point3d& queryPoint);
/**
 * @brief Actual implementation of nearestNeighbor_daniel.
 * @param[in] tree Current node in the KDTree.
 * @param[in] queryPoint Coordiantes of query point.
 * @param[in,out] minDist Pointer to currently smallest distance.
 * @param[in,out] minPoint Pointer to currently nearest point.
 * @param[in] depth Depth of the current node in the tree (needed for coordiante selection)
 */
void nearestNeighbor_daniel_impl(Node* tree, const double* queryPoint, double* minDist,
                                       Point3d* minPoint, unsigned int depth);

Point3d findNearestPoint_Elke(Node* tree, const Point3d& queryPoint);
void findNearestPointRecursiveley_Elke(Node* tree, double* currDist, unsigned int depth, Point3d* currentPoint);
void insertGivenPoint(Node* tree, double* currDist, unsigned int depth, Point3d* currentPoint);
double getDimension(Point3d point, unsigned int depth);




















/**
@brief	initializes the spherical range query
@param	tree is the K-D-Tree which stores the point cloud
@param	center is the middle of the sphere
@param	radius specifies the search area around the center point
@param	flagging is false by default and only used for efficient thinning
*/
std::vector<Point3d*> queryRangeSphere(Node* tree, Point3d* center, const double radius, bool flagging);

/**
@brief	implements the spherical range query
@param	tree is the K-D-Tree which stores the point cloud
@param	center is the middle of the sphere
@param	radius specifies the search area around the center point
@param	depth stores the current level of the traversed tree
@param	out stores the points founds in the neighborhood
@param	flagging is false by default and only used for efficient thinning
*/
void queryRangeSphere_impl(Node* tree, Point3d* center, const double radius, unsigned int depth, std::vector<Point3d*>& out, bool flagging);

/**
@brief	Thins a given pointcloud (stored in K-D-Tree) by a given radius.
@param	globalTree root node of the K-D-Tree that stores the original point data
@param	subTree node of the currently (used for recursive traversion)
@param	radius specifies the degree of thinning
@param	output stores the remaining points
*/
void homogeneousThinning(Node* globalTree, Node* subTree, const double radius, std::vector<Point3d>& output);

/** @} */

#endif // KD_TREE_H