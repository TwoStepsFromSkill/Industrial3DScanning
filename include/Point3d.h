#ifndef MY_POINT_H
#define MY_POINT_H

#include <cstdlib>

struct Point3d
{
    Point3d()
        : data{0.0, 0.0, 0.0}
    {}

    Point3d(double _x, double _y, double _z)
        : data{_x, _y, _z}
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

    double data[3];
};

double  sqr(double value);                                  ///< returns the square of a value
double  vectorLength(const Point3d& v);                     ///< returns the length of a vector
double  dotProduct  (const Point3d& v1, const Point3d& v2); ///< returns the dot product of two 3d vectors
Point3d crossProduct(const Point3d& v1, const Point3d& v2); ///< returns the cross product of two 3d vectors
void    normalizeVector(Point3d& v);                        ///< normalizes a 3d vector
double  sqDistance3d(const Point3d& v1, const Point3d& v2); ///< returns the squared Euclidean distance between two 3d points/vectors
double  distance3d  (const Point3d& v1, const Point3d& v2); ///< returns the Euclidean distance between two 3d points/vectors


#endif //end the include guard MY_POINT_H
