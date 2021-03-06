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
    Point3d operator + (const Point3d& p2) const;
    /**
     * @brief Vector subtraction.
     * @param[in] p2 Vector to subtractfrom the one stored inside the struct.
     * @return Returns a new vector that is the subtraction of the vector stored
     *          inside the struct and p2.
     */
    Point3d operator - (const Point3d& p2) const;
    /**
     * @brief Vector scalar multiplication.
     * @param[in] scalar A scalar value.
     * @returns Returns a scaled version of the vector stored inside the struct.
     */
    Point3d operator * (double scalar) const;

    /**
     * @brief Add vector onto the current one.
     * @param[in] p2 Vector to add the the one stored inside the struct.
     * @return Returns a reference to itself.
     */
    Point3d& operator += (const Point3d& p2);
    /**
     * @brief Subtract vector from the current one.
     * @param[in] p2 Vector to subtract from the one stored inside the struct.
     * @return Returns a reference to itself.
     */
    Point3d& operator -= (const Point3d& p2);
    /**
     * @brief Multiply current vector with a scalar.
     * @param[in] scalar A scalar value.
     * @returns Returns a reference to itself.
     */
    Point3d& operator *= (double scalar);

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
    bool operator == (const Point3d& p2);
    /**
     * @brief Check if vectors are inequal.
     */
    bool operator != (const Point3d& p2);

    double data[3];
    bool flag_ignore;
};

/** @defgroup group_point Point3D methods
 *  @{
 */

/**
 * @brief Square a value.
 */
double  sqr(double value);
/**
 * @brief Computes length of a given vector.
 *
 * Computes the length as \f$ s = \sqrt{p_{x}^2 + p_{y}^2 + p_{z}^2}\f$.
 */
double  vectorLength(const Point3d& v);
/**
 * @brief Computes the inner product between two vectors.
 *
 * The dot product: \f[
 *  \mathbf{u}^{T}\mathbf{v} = \sum_{i = 0}^{n}{u_{i}v_{i}}
 * \f]
 */
double  dotProduct  (const Point3d& v1, const Point3d& v2);
/**
 * @brief Computes the cross product between two vectors.
 */
Point3d crossProduct(const Point3d& v1, const Point3d& v2);
/**
 * @brief Normalizes a given vector inplace.
 * @post vectorLength(v) == 1.0
 *
 * The normalized vector is: \f[
 *  \frac{\mathbf{n}}{\|\mathbf{n}\|}
 * \f].
 */
void    normalizeVector(Point3d& v);
/**
 * @brief Computes squared distance between two points.
 *
 * This computes: \f[
 *  (v_{x} - u_{x})^2 + (v_{y} - u_{y})^2 + (v_{y} - u_{y})^2
 * \f]
 */
double  sqDistance3d(const Point3d& v1, const Point3d& v2);
/**
 * @brief Computes distance between two points.
 *
 *  * This computes: \f[
 *  \sqrt{(v_{x} - u_{x})^2 + (v_{y} - u_{y})^2 + (v_{y} - u_{y})^2}
 * \f]
 */
double  distance3d  (const Point3d& v1, const Point3d& v2);

/** @} */

#endif //end the include guard MY_POINT_H
