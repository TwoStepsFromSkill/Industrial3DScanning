#include "GLcamera.h"

#include <cmath>
#include <QtOpenGL/qgl.h>
#include <GL/glu.h>

GLcamera::GLcamera()
    : m_rotationCenter{0.0, 0.0, 0.0}
    , m_sceneRadius{1.0}
    , m_zoomFactor{1.0}
    , m_winWidth{100}
    , m_winHeight{100}
    , m_zNear{0.00001}
    , m_zFar{10.0}
    , m_usePerspectiveProjection{false}
{}

void GLcamera::initializeCamera(Point3d rotationCenter, double sceneRadius)
{
    m_rotationCenter = rotationCenter;
    m_sceneRadius = sceneRadius;
    m_zoomFactor = 1;

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Ermitteln des Abstandes der Kamera zur Szenenmitte, sodass die gesamte Szene sichtbar ist
    const double D = m_sceneRadius / sin((3.1415 / 180.0) * (45.0 / 2));

    m_zNear = D - 3 * m_sceneRadius;
    m_zFar  = D + 3 * m_sceneRadius;

    glTranslated(-rotationCenter[0], -rotationCenter[1], -rotationCenter[2] - D);
}

void GLcamera::setWindowSize(int winWidth, int winHeight)
{
    m_winWidth = winWidth;
    m_winHeight = winHeight;
}

void GLcamera::updateProjection()
{
    const double aspect = static_cast<double>(m_winWidth) / m_winHeight;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    if (m_usePerspectiveProjection)
    {
        double m_angleFOV = 45;
        double camAngle = m_angleFOV / m_zoomFactor;

        if (camAngle < 0.25)
        {
            camAngle = 0.25;
            m_zoomFactor = m_angleFOV / 0.25;
        }
        else if (camAngle > 120.0)
        {
            camAngle = 120.0;
            m_zoomFactor = m_angleFOV / 120.0;
        }

        gluPerspective(camAngle, aspect, 0.001, m_zFar);
    }
    else
    {
        const double range = m_sceneRadius / m_zoomFactor;

        if (m_winWidth <= m_winHeight)
            glOrtho(-range, range, -range / aspect, range / aspect, 0.0001, m_zFar);
        else
            glOrtho(-range*aspect, range*aspect, -range, range, 0.0001, m_zFar);
    }

    glMatrixMode(GL_MODELVIEW);
}

void GLcamera::zoom(double factor)
{
    m_zoomFactor *= factor;
    updateProjection();
}

void GLcamera::rotate(int x1, int y1, int x2, int y2)
{
    // The center of our virtual rotation ball is in the center of the screen
    const double x0 = m_winWidth / 2;
    const double y0 = m_winHeight / 2;
    // We set the radius of rotation ball to half of the screen height
    const double r = m_winHeight / 2; //ball radius is half the window height
    // const double r  = sqrt(sqr(m_winWidth) + sqr(m_winHeight)) / 2; //ball radius is half the window diagonal;

    double lastposX = x1;
    double lastposY = y1;
    double currposX = x2;
    double currposY = y2;

    double lastPosZ = (sqr(r) - sqr(lastposX - x0) - sqr(lastposY - y0));
    double currPosZ = (sqr(r) - sqr(currposX - x0) - sqr(currposY - y0));

    // if z is negative then we are outside the virtual ball and the rotation is just around the current view
    // direction/z-axis
    if (lastPosZ<0)
        lastPosZ = 0;
    else
        lastPosZ = sqrt(lastPosZ);

    if (currPosZ<0)
        currPosZ = 0;
    else
        currPosZ = sqrt(currPosZ);

    // take into account that the screen origin is in the top left corner (and not bottom left) -> x'=x-x0 and y'=y0-y
    Point3d lastPos3d(lastposX - x0, y0 - lastposY, lastPosZ);
    Point3d currPos3d(currposX - x0, y0 - currposY, currPosZ);

    normalizeVector(lastPos3d); //make unit normal vector
    normalizeVector(currPos3d); //make unit normal vector

    // the current mouse interaction results in this 3d rotation in camera space (unit sphere)
    Point3d axisCS = crossProduct(lastPos3d, currPos3d);
    double  angle = acos(dotProduct(lastPos3d, currPos3d));

    // The current rotation now needs to be multiplied with the global world rotation/transform
    // Therefore we ask OpenGL to give ous the global scene transform (which is stored in the Modelview-Matrix)
    double M[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, M); //note that OpenGL-Matrices are column-major

    // We now multiply our current rotation axis (in camera space) with the global world transform Matrix
    // and get a new rotation axis which is now in the frame of the global transform
    Point3d axisWS;
    axisWS[0] = (M[0] * axisCS[0]) + (M[1] * axisCS[1]) + (M[2] * axisCS[2]);
    axisWS[1] = (M[4] * axisCS[0]) + (M[5] * axisCS[1]) + (M[6] * axisCS[2]);
    axisWS[2] = (M[8] * axisCS[0]) + (M[9] * axisCS[1]) + (M[10] * axisCS[2]);

    // Rotation always happens in the origin (rotation center)
    // Therefore we first "move" the frame origin to our own rotation center
    glTranslated(m_rotationCenter[0], m_rotationCenter[1], m_rotationCenter[2]);
    // now we rotate the frame about our own origin/rotation center
    glRotated(2*angle * 180.0 / 3.1415, axisWS[0], axisWS[1], axisWS[2]);
    // and finally we "move" the frame origin back
    glTranslated(-m_rotationCenter[0], -m_rotationCenter[1], -m_rotationCenter[2]);
}

void GLcamera::usePerspectiveProjection(bool enable)
{
    m_usePerspectiveProjection = enable;
    updateProjection();
}

bool GLcamera::usesPerspectiveProjection() const
{
    return m_usePerspectiveProjection;
}
