//the following #-syntax at the begin AND end of this file is an "include guard". 
//It ensures that the compiler include/uses this file only once (if it used multiple times)
#ifndef MY_POINT_H
#define MY_POINT_H

//Definition of Point3d structure (our 3d Vector)
struct Point3d
{
public:
    Point3d() /// < constructor of this object
    {
        x = 0; y = 0; z = 0;
    }
    Point3d(double _x, double _y, double _z) ///< another constructor of this class
    {
        x = _x; y = _y; z = _z;
    }

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

    double x,y,z;
};

double  sqr(double value);                                  ///< returns the square of a value
double  vectorLength(const Point3d& v);                     ///< returns the length of a vector
double  dotProduct  (const Point3d& v1, const Point3d& v2); ///< returns the dot product of two 3d vectors
Point3d crossProduct(const Point3d& v1, const Point3d& v2); ///< returns the cross product of two 3d vectors
void    normalizeVector(Point3d& v);                        ///< normalizes a 3d vector
double  sqDistance3d(const Point3d& v1, const Point3d& v2); ///< returns the squared Euclidean distance between two 3d points/vectors
double  distance3d  (const Point3d& v1, const Point3d& v2); ///< returns the Euclidean distance between two 3d points/vectors


#endif //end the include guard MY_POINT_H
