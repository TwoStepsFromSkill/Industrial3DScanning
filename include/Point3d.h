#ifndef MY_POINT_H
#define MY_POINT_H

#include <cstdlib>

/**
 * @brief Point in 3-dimensional euclidean space.
 */
struct Point3d
{
    /**
     * @brief Default constructor. Initializes to 0.
     */
    Point3d()
        : data{0.0, 0.0, 0.0}, flag_ignore{false}
    {}

    /**
     * @brief Constructor. Initializes with given x,y,z.
     * @param x,y,z Coordinates of point.
     */
    Point3d(double _x, double _y, double _z)
        : data{_x, _y, _z}, flag_ignore{false}
    {}

    //custom operators that enable vector algebra
    //these operator are marked CONST because they don't change member variables of this class
    Point3d operator + (const Point3d& p2) const; // + operation of two points (pt3 = pt1 + pt2)
    Point3d operator - (const Point3d& p2) const; // - operation of two points (pt3 = pt1 - pt2)
    Point3d operator * (double scalar) const;     // multiplication with a scalar value(pt2 = pt1 * 0.5)

    //assignments with operator
    //the operators can not be marked CONST because they do change the input
    Point3d& operator += (const Point3d& p2); // += operation of two points (pt1+= pt2  -> pt1 = pt1 + pt2)
    Point3d& operator -= (const Point3d& p2); // += operation of two points (pt1-= pt2  -> pt1 = pt1 - pt2)
    Point3d& operator *= (double scalar);     // *= multiplication with a scalar (pt1*= s  -> pt1 = pt1 * s)

    double& operator [] (std::size_t index);
    const double& operator [] (std::size_t index) const;

    bool operator == (const Point3d& p2);	// == comparison of two points (p1 == p2 -> true)
    bool operator != (const Point3d& p2);	// != comparison of two points (p1 != p2 -> true)

    double data[3];
    bool flag_ignore;
};

double  sqr(double value);                                  ///< returns the square of a value
double  vectorLength(const Point3d& v);                     ///< returns the length of a vector
double  dotProduct  (const Point3d& v1, const Point3d& v2); ///< returns the dot product of two 3d vectors
Point3d crossProduct(const Point3d& v1, const Point3d& v2); ///< returns the cross product of two 3d vectors
void    normalizeVector(Point3d& v);                        ///< normalizes a 3d vector
double  sqDistance3d(const Point3d& v1, const Point3d& v2); ///< returns the squared Euclidean distance between two 3d points/vectors
double  distance3d  (const Point3d& v1, const Point3d& v2); ///< returns the Euclidean distance between two 3d points/vectors


#endif //end the include guard MY_POINT_H
