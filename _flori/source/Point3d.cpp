#include <cmath>        //for standard C/C++ math functions
#include "Point3d.h" 

//custom operator that enables the + operation of two points (pt3 = pt1 + pt2)
Point3d Point3d::operator + (const Point3d& p2) const
{
  Point3d result;
  result.x = x + p2.x;
  result.y = y + p2.y;
  result.z = z + p2.z;
  return result;
}

//custom operator that enables the - operation of two points (pt3 = pt1 - pt2)
Point3d Point3d::operator - (const Point3d& p2) const
{
  Point3d result;
  result.x = x - p2.x;
  result.y = y - p2.y;
  result.z = z - p2.z;
  return result;
}
//custom operator that enables the multiplication with a scalar value (pt2 = pt1 * 0.5)
Point3d Point3d::operator * (double scalar) const
{
  Point3d result;
  result.x = scalar * x;
  result.y = scalar * y;
  result.z = scalar * z;
  return result;
}

//custom operator that enables the += operation (pt1 += pt2 -> pt1 = pt1 + pt2)
Point3d& Point3d::operator += (const Point3d& p2)
{
  x += p2.x;
  y += p2.y;
  z += p2.z;
  return *this;
}

//custom operator that enables the -= operation (pt1 -= pt2 -> pt1 = pt1 - pt2)
Point3d& Point3d::operator -= (const Point3d& p2)
{
  x -= p2.x;
  y -= p2.y;
  z -= p2.z;
  return *this;
}

//custom operator that enables the += operation (pt1 *= 2 -> pt1 = pt1 * s)
Point3d& Point3d::operator *= (double scalar)
{
  x *= scalar;
  y *= scalar;
  z *= scalar;
  return *this;
}

//custom operator that enables the == comparison between two points
bool Point3d::operator == (Point3d& p2)
{
	return x == p2.x && y == p2.y && z == p2.z;
}

//custom operator that enables the != comparison between two points
bool Point3d::operator != (Point3d& p2)
{
	return x != p2.x | y != p2.y | z != p2.z;
}















//returns the square of a value (unfortunately C++ does not provide this function itself...)
double sqr(double value)
{
  return value*value;
}

//returns the length of a vector
double  vectorLength(const Point3d& v)
{
  double length = sqrt(sqr(v.x) + sqr(v.y) + sqr(v.z));
  return length;
}

//returns the dot product of two 3d vectors
double dotProduct(const Point3d& v1, const Point3d& v2)
{
  return (v1.x*v2.x) + (v1.y*v2.y) + (v1.z*v2.z);
}

//returns the cross product of two 3d vectors
Point3d crossProduct(const Point3d& v1, const Point3d& v2)
{
  Point3d result;
  result.x = (v1.y * v2.z) - (v1.z * v2.y);
  result.y = (v1.z * v2.x) - (v1.x * v2.z);
  result.z = (v1.x * v2.y) - (v1.y * v2.x);

  return result;
}

//normalizes a 3d vector (direction vector)
void normalizeVector(Point3d& v)
{
  const double length = vectorLength(v);
  if (length > 0)
  {
    v.x /= length;
    v.y /= length;
    v.z /= length;
  }
}

///< returns the squared Euclidean distance between two 3d points/vectors
double sqDistance3d(const Point3d& v1, const Point3d& v2)
{
  const double d = sqr(v1.x - v2.x) + sqr(v1.y - v2.y) + sqr(v1.z - v2.z);
  return d;
}

//returns the Euclidean distance between two 3d points / vectors
double distance3d(const Point3d& v1, const Point3d& v2)
{
  const double d = std::sqrt(sqDistance3d(v1,v2));
  return d;
}

//returns the point coordinates as string
std::string Point3dToString(Point3d* pt)
{
	return "(" + std::to_string(pt->x) + ", " + std::to_string(pt->y) + ", " + std::to_string(pt->z) + ")";
}
