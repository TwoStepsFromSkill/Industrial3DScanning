#include <cmath>        //for standard C/C++ math functions
#include "Point3d.h"

double& Point3d::operator [] (std::size_t index)
{
    return data[index];
}

const double& Point3d::operator [] (std::size_t index) const
{
    return data[index];
}

//custom operator that enables the + operation of two points (pt3 = pt1 + pt2)
Point3d Point3d::operator + (const Point3d& p2) const
{
    return Point3d(data[0] + p2[0], data[1] + p2[1], data[2] + p2[2]);
}

//custom operator that enables the - operation of two points (pt3 = pt1 - pt2)
Point3d Point3d::operator - (const Point3d& p2) const
{
    return Point3d(data[0] - p2[0], data[1] - p2[1], data[2] - p2[2]);
}
//custom operator that enables the multiplication with a scalar value (pt2 = pt1 * 0.5)
Point3d Point3d::operator * (double scalar) const
{
    return Point3d(data[0] * scalar, data[1] * scalar, data[2] * scalar);
}

//custom operator that enables the += operation (pt1 += pt2 -> pt1 = pt1 + pt2)
Point3d& Point3d::operator += (const Point3d& p2)
{
    data[0] += p2[0];
    data[1] += p2[1];
    data[2] += p2[1];
    return *this;
}

//custom operator that enables the -= operation (pt1 -= pt2 -> pt1 = pt1 - pt2)
Point3d& Point3d::operator -= (const Point3d& p2)
{
    data[0] -= p2[0];
    data[1] -= p2[1];
    data[2] -= p2[1];
    return *this;
}

//custom operator that enables the += operation (pt1 *= 2 -> pt1 = pt1 * s)
Point3d& Point3d::operator *= (double scalar)
{
    data[0] *= scalar;
    data[1] *= scalar;
    data[2] *= scalar;
    return *this;
}

//returns the square of a value (unfortunately C++ does not provide this function itself...)
double sqr(double value)
{
    return value*value;
}

//returns the length of a vector
double  vectorLength(const Point3d& v)
{
    double length = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
    return length;
}

//returns the dot product of two 3d vectors
double dotProduct(const Point3d& v1, const Point3d& v2)
{
    return (v1[0]*v2[0]) + (v1[1]*v2[1]) + (v1[2]*v2[2]);
}

//returns the cross product of two 3d vectors
Point3d crossProduct(const Point3d& v1, const Point3d& v2)
{
    return Point3d((v1[1] * v2[2]) - (v1[2] * v2[1]),
                   (v1[2] * v2[0]) - (v1[0] * v2[2]),
                   (v1[0] * v2[1]) - (v1[1] * v2[0]));
}

//normalizes a 3d vector (direction vector)
void normalizeVector(Point3d& v)
{
    const double length = vectorLength(v);
    if (length > 0)
    {
        v[0] /= length;
        v[1] /= length;
        v[2] /= length;
    }
}

///< returns the squared Euclidean distance between two 3d points/vectors
double sqDistance3d(const Point3d& v1, const Point3d& v2)
{
    return sqr(v1[0] - v2[0]) + sqr(v1[1] - v2[1]) + sqr(v1[2] - v2[2]);
}

//returns the Euclidean distance between two 3d points / vectors
double distance3d(const Point3d& v1, const Point3d& v2)
{
    return std::sqrt(sqDistance3d(v1,v2));
}
