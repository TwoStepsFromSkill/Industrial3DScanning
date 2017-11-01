#include "GLwidget.h"
#include <GL/glu.h>
#include <QtGui/QMouseEvent>
#include <QSettings>
#include <QVariant>
#include <QFileInfo>

#include <cmath>

GLwidget::GLwidget(QWidget* parent)
    : QOpenGLWidget(parent)
    , m_points()
    , m_pointsInRange()
    , m_mouseLastPos()
    , m_camera()
    , m_bbmin()
    , m_bbmax()
    , m_sceneCenter()
    , m_sceneRadius()
    , m_rangeCenter()
    , m_rangeExtend()
    , m_drawRangeQueryBox(false)
    , m_drawRangeQueryResult(false)
{
    loadDrawSettings();
}

void GLwidget::initializeGL()
{
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void GLwidget::resizeGL(int w, int h)
{
    glViewport(0,0,w,h);

    m_camera.setWindowSize(w,h);
    m_camera.updateProjection(); //adjust projection to new window size
}

void GLwidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); //clear buffers
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);   //clear background color
    glClearDepth(1.0f);                     //clear depth buffer

    //draws the scene background
    drawBackground();

    //Draw pointclouds
    if (!m_points.empty())
    { /* Drawing Points with VertexArrays */
        glEnableClientState(GL_VERTEX_ARRAY);

        glPointSize(m_PC_size);
        glColor3ubv(m_PC_color);
        glVertexPointer(3, GL_DOUBLE, sizeof(Point3d), &m_points[0]);
        glDrawArrays(GL_POINTS, 0, (unsigned int)m_points.size());

        glDisableClientState(GL_VERTEX_ARRAY);
    }


    if (m_drawRangeQueryBox)
    {
        // Draw range box
        glPushAttrib(GL_POLYGON_BIT);
        glColor3ubv(m_RQ_box_color);
        glLineWidth(m_RQ_box_width);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

        glBegin(GL_QUADS);
            glVertex3d(m_fbl[0], m_fbl[1], m_fbl[2]); glVertex3d(m_fbr[0], m_fbr[1], m_fbr[2]);
            glVertex3d(m_ftr[0], m_ftr[1], m_ftr[2]); glVertex3d(m_ftl[0], m_ftl[1], m_ftl[2]);  // front

            glVertex3d(m_bbl[0], m_bbl[1], m_bbl[2]); glVertex3d(m_bbr[0], m_bbr[1], m_bbr[2]);
            glVertex3d(m_btr[0], m_btr[1], m_btr[2]); glVertex3d(m_btl[0], m_btl[1], m_btl[2]);  // back

            glVertex3d(m_fbl[0], m_fbl[1], m_fbl[2]); glVertex3d(m_fbr[0], m_fbr[1], m_fbr[2]);
            glVertex3d(m_bbr[0], m_bbr[1], m_bbr[2]); glVertex3d(m_bbl[0], m_bbl[1], m_bbl[2]);  // bottom

            glVertex3d(m_ftl[0], m_ftl[1], m_ftl[2]); glVertex3d(m_ftr[0], m_ftr[1], m_ftr[2]);
            glVertex3d(m_btr[0], m_btr[1], m_btr[2]); glVertex3d(m_btl[0], m_btl[1], m_btl[2]);  // top

            glVertex3d(m_ftl[0], m_ftl[1], m_ftl[2]); glVertex3d(m_fbl[0], m_fbl[1], m_fbl[2]);
            glVertex3d(m_bbl[0], m_bbl[1], m_bbl[2]); glVertex3d(m_btl[0], m_btl[1], m_btl[2]);  // left

            glVertex3d(m_ftr[0], m_ftr[1], m_ftr[2]); glVertex3d(m_fbr[0], m_fbr[1], m_fbr[2]);
            glVertex3d(m_bbr[0], m_bbr[1], m_bbr[2]); glVertex3d(m_btr[0], m_btr[1], m_btr[2]);  // right

        glEnd();
        glLineWidth(1.0f);
        glPopAttrib();
    }

    if (m_drawRangeQueryResult)
    {
        // Draw points
        glEnableClientState(GL_VERTEX_ARRAY);

        glPointSize(m_RQ_result_size);
        glColor3ubv(m_RQ_result_color);
        glVertexPointer(3, GL_DOUBLE, sizeof(Point3d), &m_pointsInRange[0]);
        glDrawArrays(GL_POINTS, 0, (unsigned int)m_pointsInRange.size());

        glPointSize(m_PC_size);
        glDisableClientState(GL_VERTEX_ARRAY);
    }

    glEnable(GL_DEPTH_TEST);

    //draw coordinate frame
    drawCoordinateAxes();

    update();
}

void  GLwidget::setPointsInRange(const std::vector<Point3d>& points)
{
    m_pointsInRange = points;
    this->repaint();
}

void GLwidget::drawingRangeQueryBoxChange(bool value)
{
    m_drawRangeQueryBox = value;
    m_drawRangeQueryResult = value;

    this->repaint();
}

void GLwidget::drawingRangeQueryResultEnabled(bool value)
{
    m_drawRangeQueryResult = value;

    if (!value)
        m_pointsInRange.clear();

    this->repaint();
}

void GLwidget::rangeQueryCenterChanged(double x, double y, double z)
{
    m_rangeCenter.x = x;
    m_rangeCenter.y = y;
    m_rangeCenter.z = z;

    updateRangeQueryBoxData();
    this->repaint();
}

void GLwidget::rangeQueryExtendChanged(double dx, double dy, double dz)
{
    m_rangeExtend.x = dx;
    m_rangeExtend.y = dy;
    m_rangeExtend.z = dz;

    updateRangeQueryBoxData();
    this->repaint();
}

void GLwidget::loadDrawSettings()
{
    QFileInfo check_file(QString("drawSettings.ini"));

    if (!(check_file.exists() && check_file.isFile()))
        std::cerr << "Could not find drawSettings.ini file in executable directory. Use default values!\n";

    QSettings set(QString("drawSettings.ini"), QSettings::IniFormat);

    m_PC_color[0] = set.value("PC_color_R", 255).value<unsigned char>();
    m_PC_color[1] = set.value("PC_color_G", 141).value<unsigned char>();
    m_PC_color[2] = set.value("PC_color_B", 42).value<unsigned char>();
    m_PC_size = set.value("PC_size", 2).toInt();

    m_RQ_box_color[0] = set.value("RQ_box_color_R", 188).value<unsigned char>();
    m_RQ_box_color[1] = set.value("RQ_box_color_G", 217).value<unsigned char>();
    m_RQ_box_color[2] = set.value("RQ_box_color_B", 5).value<unsigned char>();
    m_RQ_box_width = set.value("RQ_box_width", 1.0).toDouble();

    m_RQ_result_color[0] = set.value("RQ_result_color_R", 217).value<unsigned char>();
    m_RQ_result_color[1] = set.value("RQ_result_color_G", 22).value<unsigned char>();
    m_RQ_result_color[2] = set.value("RQ_result_color_B", 25).value<unsigned char>();
    m_RQ_result_size = set.value("RQ_result_size", 2).toInt();

    m_XAXIS_color[0] = set.value("XAXIS_color_R", 255).value<unsigned char>();
    m_XAXIS_color[1] = set.value("XAXIS_color_G", 0).value<unsigned char>();
    m_XAXIS_color[2] = set.value("XAXIS_color_B", 0).value<unsigned char>();

    m_YAXIS_color[0] = set.value("YAXIS_color_R", 0).value<unsigned char>();
    m_YAXIS_color[1] = set.value("YAXIS_color_G", 255).value<unsigned char>();
    m_YAXIS_color[2] = set.value("YAXIS_color_B", 0).value<unsigned char>();

    m_ZAXIS_color[0] = set.value("ZAXIS_color_R", 0).value<unsigned char>();
    m_ZAXIS_color[1] = set.value("ZAXIS_color_G", 0).value<unsigned char>();
    m_ZAXIS_color[2] = set.value("ZAXIS_color_B", 255).value<unsigned char>();
    m_AXIS_width = set.value("AXIS_width", 1.0).toDouble();

    m_CENTERSPHERE_color[0] = set.value("CENTERSPHERE_color_R", 255).value<unsigned char>();
    m_CENTERSPHERE_color[1] = set.value("CENTERSPHERE_color_G", 255).value<unsigned char>();
    m_CENTERSPHERE_color[2] = set.value("CENTERSPHERE_color_B", 0).value<unsigned char>();

    m_XCIRCLE_color[0] = set.value("XCIRCLE_color_R", 255).value<unsigned char>();
    m_XCIRCLE_color[1] = set.value("XCIRCLE_color_G", 0).value<unsigned char>();
    m_XCIRCLE_color[2] = set.value("XCIRCLE_color_B", 0).value<unsigned char>();

    m_YCIRCLE_color[0] = set.value("YCIRCLE_color_R", 0).value<unsigned char>();
    m_YCIRCLE_color[1] = set.value("YCIRCLE_color_G", 255).value<unsigned char>();
    m_YCIRCLE_color[2] = set.value("YCIRCLE_color_B", 0).value<unsigned char>();

    m_ZCIRCLE_color[0] = set.value("ZCIRCLE_color_R", 0).value<unsigned char>();
    m_ZCIRCLE_color[1] = set.value("ZCIRCLE_color_G", 0).value<unsigned char>();
    m_ZCIRCLE_color[2] = set.value("ZCIRCLE_color_B", 255).value<unsigned char>();
    m_CIRCLE_width = set.value("CIRCLE_width", 1.0).toDouble();

    m_BB_color[0] = set.value("BB_color_R", 255).value<unsigned char>();
    m_BB_color[1] = set.value("BB_color_G", 255).value<unsigned char>();
    m_BB_color[2] = set.value("BB_color_B", 255).value<unsigned char>();
    m_BB_width = set.value("BB_width", 1.0).toDouble();

    m_BG_top_color[0] = set.value("BG_top_color_R", 100).value<unsigned char>();
    m_BG_top_color[1] = set.value("BG_top_color_G", 100).value<unsigned char>();
    m_BG_top_color[2] = set.value("BG_top_color_B", 100).value<unsigned char>();

    m_BG_bottom_color[0] = set.value("BG_bottom_color_R", 38).value<unsigned char>();
    m_BG_bottom_color[1] = set.value("BG_bottom_color_G", 38).value<unsigned char>();
    m_BG_bottom_color[2] = set.value("BG_bottom_color_B", 38).value<unsigned char>();
}

void GLwidget::mousePressEvent(QMouseEvent * e)  ///<
{
    if (e->buttons() == Qt::LeftButton)
        m_mouseLastPos = e->pos();
}

void GLwidget::mouseMoveEvent(QMouseEvent * e)   ///<
{
    if (e->buttons() != Qt::LeftButton){ return; }

    //No action, if mouse position outside the window
    if ((*e).x()<0 || (*e).x() >= width() ){ return; }
    if ((*e).y()<0 || (*e).y() >= height()){ return; }

    const QPoint& mouseCurrentPos(e->pos());

    if (m_mouseLastPos == mouseCurrentPos) return;

    makeCurrent();
    m_camera.rotate(m_mouseLastPos.x(), m_mouseLastPos.y(), mouseCurrentPos.x(), mouseCurrentPos.y());

    //the current mouse position is the last mouse position for the next interaction
    m_mouseLastPos = mouseCurrentPos;

    update();
}

void GLwidget::wheelEvent(QWheelEvent * e)
{
    if (e->orientation() != Qt::Vertical) return;

    const double factor = e->delta()<0 ? 1.1 : 0.9;

    makeCurrent();
    m_camera.zoom(factor);
}

void GLwidget::updateScene()
{
    if (m_points.empty())
    {
        return;
    }

    //OK, we now compute the min and max coordinates for our bounding box
    Point3d minPoint = m_points.front(); //initialize min with the first point
    Point3d maxPoint = m_points.front(); //initialize max with the first point

    for (unsigned int i = 0; i < m_points.size(); ++i)
    {
        const Point3d& pt = m_points[i]; //do not copy but get a reference to the i-th point in the vector
        if (pt.x < minPoint.x) minPoint.x = pt.x;
        else if (pt.x > maxPoint.x) maxPoint.x = pt.x;

        if (pt.y < minPoint.y) minPoint.y = pt.y;
        else if (pt.y > maxPoint.y) maxPoint.y = pt.y;

        if (pt.z < minPoint.z) minPoint.z = pt.z;
        else if (pt.z > maxPoint.z) maxPoint.z = pt.z;
    }

    m_sceneCenter.x = (maxPoint.x + minPoint.x) / 2;
    m_sceneCenter.y = (maxPoint.y + minPoint.y) / 2;
    m_sceneCenter.z = (maxPoint.z + minPoint.z) / 2;

    m_sceneRadius = distance3d(m_sceneCenter, maxPoint);

    m_bbmin=minPoint;
    m_bbmax=maxPoint;

//     std::cout << "\nBounding Box was computed:\n";
//     std::cout << "minPoint is: " << minPoint.x << "," << minPoint.y << "," << minPoint.z << std::endl;
//     std::cout << "maxPoint is: " << maxPoint.x << "," << maxPoint.y << "," << maxPoint.z << std::endl;

    makeCurrent();
    m_camera.initializeCamera(m_sceneCenter, m_sceneRadius); //set rotation center and scene radius to initially setup the camera
    m_camera.updateProjection(); //updateProjection to make the new point cloud fit to the screen
}

/** draws a unit box.
origin is at (0,0,0) and maximum at (1,1,1).
use glTranslate to move the origin to the minimum coordinates
afterwards use glScale(sx,sy,sz) resize the box.
*/
void GLwidget::drawBox()
{
    glBegin(GL_QUADS);
    glVertex3d(0, 0, 0); glVertex3d(1, 0, 0); glVertex3d(1, 1, 0); glVertex3d(0, 1, 0); //front
    glVertex3d(0, 0, 1); glVertex3d(1, 0, 1); glVertex3d(1, 1, 1); glVertex3d(0, 1, 1); //back
    glVertex3d(0, 0, 0); glVertex3d(0, 0, 1); glVertex3d(0, 1, 1); glVertex3d(0, 1, 0); //left
    glVertex3d(1, 0, 0); glVertex3d(1, 0, 1); glVertex3d(1, 1, 1); glVertex3d(1, 1, 0); //right
    glVertex3d(0, 0, 0); glVertex3d(1, 0, 0); glVertex3d(1, 0, 1); glVertex3d(0, 0, 1); //bottom
    glVertex3d(0, 1, 0); glVertex3d(1, 1, 0); glVertex3d(1, 1, 1); glVertex3d(0, 1, 1); //top
    glEnd();
}

//draws a unit circle with radius 1 at the origin(0,0,0)
//use glTranslate, glScale and glRotate to resize and reposition the circle
void GLwidget::drawCircle()
{
    const int segments = 180;

    glBegin(GL_LINE_LOOP);
    glLineWidth(m_CIRCLE_width);
    for (int i = 0; i < segments; ++i)
    {
        const double theta = 2.0 * 3.1415926 * double(i) / double(segments);
        glVertex2d(cos(theta), sin(theta));
    }
    glEnd();
}

//draws the coordinate system
void GLwidget::drawCoordinateAxes()
{
    //draw coordinate frame
    glBegin(GL_LINES);
    glLineWidth(m_AXIS_width);
    //draw line for X-Axis
    glColor3ubv(m_XAXIS_color);
    glVertex3d(m_sceneCenter.x, m_sceneCenter.y, m_sceneCenter.z);
    glVertex3d(m_sceneCenter.x + m_sceneRadius, m_sceneCenter.y, m_sceneCenter.z);
    //draw line for Y-Axis
    glColor3ubv(m_YAXIS_color);
    glVertex3d(m_sceneCenter.x, m_sceneCenter.y, m_sceneCenter.z);
    glVertex3d(m_sceneCenter.x, m_sceneCenter.y + m_sceneRadius, m_sceneCenter.z);
    //draw line for Z-Axis
    glColor3ubv(m_ZAXIS_color);
    glVertex3d(m_sceneCenter.x, m_sceneCenter.y, m_sceneCenter.z);
    glVertex3d(m_sceneCenter.x, m_sceneCenter.y, m_sceneCenter.z + m_sceneRadius);
    glEnd();

    //draw center point as a sphere
    glPushMatrix();
    glColor3ubv(m_CENTERSPHERE_color);
    glTranslated(m_sceneCenter.x, m_sceneCenter.y, m_sceneCenter.z);
    GLUquadric* quad = gluNewQuadric();
    gluSphere(quad, m_sceneRadius / 20, 30, 30);
    gluDeleteQuadric(quad);
    glPopMatrix();

    glPushMatrix();
    glColor3ubv(m_ZCIRCLE_color);
    glTranslated(m_sceneCenter.x, m_sceneCenter.y, m_sceneCenter.z);
    glScaled(m_sceneRadius, m_sceneRadius, m_sceneRadius);
    drawCircle();
    //draw another circle 90 degree rotated
    glRotated(90, 1, 0, 0);
    glColor3ubv(m_YCIRCLE_color);
    drawCircle();
    glRotated(90, 0, 1, 0);
    glColor3ubv(m_XCIRCLE_color);
    drawCircle();
    glPopMatrix();

    //draw a disk to indicate the scene radius
    /*glPushMatrix();
  glColor3ub(255, 0, 255);
  glTranslated(m_sceneCenter.x, m_sceneCenter.y, m_sceneCenter.z);
  GLUquadric* quad2 = gluNewQuadric();
  gluDisk(quad2, m_sceneRadius, m_sceneRadius*1.001, 100, 100);
  gluDeleteQuadric(quad2);*/

    //draw bounding box
    glPushMatrix();
    glPushAttrib(GL_POLYGON_BIT);
    glColor3ubv(m_BB_color);
    Point3d S=m_bbmax-m_bbmin;
    glTranslated(m_bbmin.x, m_bbmin.y, m_bbmin.z);
    glScaled(S.x,S.y,S.z);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); //draw wire frame instead of filled quads
    glLineWidth(m_BB_width);
    drawBox();
    glPopAttrib();
    glPopMatrix();
}

//------------------------------------------------------------------------------
/** draws a static permanent 2d color gradient.
    Drawing static background means that we draw a 2D rectangle, which is not
    influenced by scene rotation and by the camera projection.
    Basically we just want to draw a 2D texture/image.
*/
//------------------------------------------------------------------------------
void GLwidget::drawBackground()
{
    makeCurrent(); //each time when we call a glXXX-Funktion QT needs the current context

    const float winWidth((float)width());
    const float winHeight((float)height());

    glPushAttrib(GL_DEPTH_BUFFER_BIT | GL_LIGHTING_BIT | GL_CURRENT_BIT);

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);

    glMatrixMode(GL_PROJECTION);           //At first we select the Projection matrix
    glPushMatrix();                        //save Projektion matrix to restore it at the end
    glLoadIdentity();                      //initialize projection matrix
    gluOrtho2D(0, winWidth, 0, winHeight); //select orthographic projecton

    glMatrixMode(GL_MODELVIEW);            //now change to Modelview because we want to draw something
    glPushMatrix();
    glLoadIdentity();

    glBegin(GL_QUADS);
    glColor3ubv(m_BG_bottom_color); //color bottom
    glVertex2f(0.0f, 0.0f);  glVertex2f(winWidth, 0.0f);
    glColor3ubv(m_BG_top_color);  //color top
    glVertex2f(winWidth, winHeight);  glVertex2f(0.0f, winHeight);
    glEnd();

    glMatrixMode(GL_PROJECTION); //select to Projektionmatrix
    glPopMatrix();  //reset 2D-projection

    glMatrixMode(GL_MODELVIEW); //select MODELVIEW
    glPopMatrix();  //reset ModelviewMatrix to last state

    glPopAttrib();  //restore last attributes
}

void GLwidget::updateRangeQueryBoxData()
{
    // Front

    m_fbl[0] = m_rangeCenter.x - (m_rangeExtend.x / 2.0);
    m_fbl[1] = m_rangeCenter.y - (m_rangeExtend.y / 2.0);
    m_fbl[2] = m_rangeCenter.z - (m_rangeExtend.z / 2.0);

    m_fbr[0] = m_rangeCenter.x + (m_rangeExtend.x / 2.0);
    m_fbr[1] = m_rangeCenter.y - (m_rangeExtend.y / 2.0);
    m_fbr[2] = m_rangeCenter.z - (m_rangeExtend.z / 2.0);

    m_ftl[0] = m_rangeCenter.x - (m_rangeExtend.x / 2.0);
    m_ftl[1] = m_rangeCenter.y - (m_rangeExtend.y / 2.0);
    m_ftl[2] = m_rangeCenter.z + (m_rangeExtend.z / 2.0);

    m_ftr[0] = m_rangeCenter.x + (m_rangeExtend.x / 2.0);
    m_ftr[1] = m_rangeCenter.y - (m_rangeExtend.y / 2.0);
    m_ftr[2] = m_rangeCenter.z + (m_rangeExtend.z / 2.0);

    // Back

    m_bbl[0] = m_rangeCenter.x - (m_rangeExtend.x / 2.0);
    m_bbl[1] = m_rangeCenter.y + (m_rangeExtend.y / 2.0);
    m_bbl[2] = m_rangeCenter.z - (m_rangeExtend.z / 2.0);

    m_bbr[0] = m_rangeCenter.x + (m_rangeExtend.x / 2.0);
    m_bbr[1] = m_rangeCenter.y + (m_rangeExtend.y / 2.0);
    m_bbr[2] = m_rangeCenter.z - (m_rangeExtend.z / 2.0);

    m_btl[0] = m_rangeCenter.x - (m_rangeExtend.x / 2.0);
    m_btl[1] = m_rangeCenter.y + (m_rangeExtend.y / 2.0);
    m_btl[2] = m_rangeCenter.z + (m_rangeExtend.z / 2.0);

    m_btr[0] = m_rangeCenter.x + (m_rangeExtend.x / 2.0);
    m_btr[1] = m_rangeCenter.y + (m_rangeExtend.y / 2.0);
    m_btr[2] = m_rangeCenter.z + (m_rangeExtend.z / 2.0);
}