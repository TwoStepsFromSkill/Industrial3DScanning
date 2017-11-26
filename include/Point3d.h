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

    /**
     * @brief Vector addition.
     * @param[in] p2 Vector to add the the one stored inside the struct.
     * @return Returns a new vector that is the addition of the vector stored
     *          inside the struct and p2.
     */
    Point3d operator + (const Point3d& p2) const; // + operation of two points (pt3 = pt1 + pt2)
    /**
     * @brief Vector subtraction.
     * @param[in] p2 Vector to subtractfrom the one stored inside the struct.
     * @return Returns a new vector that is the subtraction of the vector stored
     *          inside the struct and p2.
     */
    Point3d operator - (const Point3d& p2) const; // - operation of two points (pt3 = pt1 - pt2)
    /**
     * @brief Vector scalar multiplication.
     * @param[in] scalar A scalar value.
     * @returns Returns a scaled version of the vector stored inside the struct.
     */
    Point3d operator * (double scalar) const;     // multiplication with a scalar value(pt2 = pt1 * 0.5)

    /**
     * @brief Add vector onto the current one.
     * @param[in] p2 Vector to add the the one stored inside the struct.
     * @return Returns a reference to itself.
     */
    Point3d& operator += (const Point3d& p2); // += operation of two points (pt1+= pt2  -> pt1 = pt1 + pt2)
    /**
     * @brief Subtract vector from the current one.
     * @param[in] p2 Vector to subtract from the one stored inside the struct.
     * @return Returns a reference to itself.
     */
    Point3d& operator -= (const Point3d& p2); // += operation of two points (pt1-= pt2  -> pt1 = pt1 - pt2)
    /**
     * @brief Multiply current vector with a scalar.
     * @param[in] scalar A scalar value.
     * @returns Returns a reference to itself.
     */
    Point3d& operator *= (double scalar);     // *= multiplication with a scalar (pt1*= s  -> pt1 = pt1 * s)

    /**
     * @brief Element access.
     * @param[in] index Index of the coordinate.
     * @note This method does not validate input!
     */
    double& operator [] (std::size_t index);
    /**
     * @copydoc Point3d::operator[]
     */
    const double& operator [] (std::size_t index) const;

    /**
     * @brief Compare vector to another.
     * @return Returns true iff all elements are equal.
     */
    bool operator == (const Point3d& p2);	// == comparison of two points (p1 == p2 -> true)
    /**
     * @brief Check if vectors are inequal.
     */
    bool operator != (const Point3d& p2);	// != comparison of two points (p1 != p2 -> true)

    double data[3];
    bool flag_ignore;
};

/**
 * @brief Square a value.
 */
double  sqr(double value);
/**
 * @brief Computes length of a given vector.
 */
double  vectorLength(const Point3d& v);
/**
 * @brief Computes the inner product between two vectors.
 */
double  dotProduct  (const Point3d& v1, const Point3d& v2);
/**
 * @brief Computes the cross product between two vectors.
 */
Point3d crossProduct(const Point3d& v1, const Point3d& v2);
/**
 * @brief Normalizes a given vector inplace.
 * @post vectorLength(v) == 1.0
 */
void    normalizeVector(Point3d& v);
/**
 * @brief Computes squared distance between two points.
 */
double  sqDistance3d(const Point3d& v1, const Point3d& v2);
/**
 * @brief Computes distance between two points.
 */
double  distance3d  (const Point3d& v1, const Point3d& v2);

#endif //end the include guard MY_POINT_H
