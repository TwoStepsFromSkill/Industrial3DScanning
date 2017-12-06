#include "GLwidget.h"
#include <GL/glu.h>
#include <QtGui/QMouseEvent>
#include <QSettings>
#include <QVariant>
#include <QFileInfo>
#include <QEvent>
#include <QKeyEvent>
#include <array>
#include <vector>

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
    , m_drawBFP(false)
    , m_drawBFL(false)
    , m_drawRangeQueryBox(false)
    , m_drawRangeQueryResult(false)
    , m_drawMainPointCloud(true)
    , m_drawMainPointCloudWithColorArray(false)
    , m_drawSmoothedPointCloudWithColorArray(false)
    , m_drawSmoothedState(0)
    , m_drawTemporary(false)
    , m_drawNearestQueryPoint(false)
    , m_drawNearestResultPoint(false)
    , m_currentModeIsSmoothed(false)
    , m_drawSmoothedPoints(false)
    , m_drawThinnedPoints(false)
{
    loadDrawSettings();
}

void GLwidget::initializeGL()
{
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glEnable(GL_DEPTH_TEST);

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
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
    if (!m_points.empty() && m_drawMainPointCloud)
    {
        glEnableClientState(GL_VERTEX_ARRAY);

        glPointSize(m_PC_size);

        if (m_drawMainPointCloudWithColorArray)
        {
            glEnableClientState(GL_COLOR_ARRAY);
            glVertexPointer(3, GL_DOUBLE, sizeof(Point3d), &m_points[0]);
            glColorPointer(3, GL_UNSIGNED_BYTE, 3*sizeof(unsigned char), &m_pointColors[0]);
            glDrawArrays(GL_POINTS, 0, (unsigned int) m_points.size());

            glDisableClientState(GL_COLOR_ARRAY);
        }
        else
        {
            glColor4ubv(m_PC_color);
            glVertexPointer(3, GL_DOUBLE, sizeof(Point3d), &m_points[0]);
            glDrawArrays(GL_POINTS, 0, (unsigned int) m_points.size());
        }

        glDisableClientState(GL_VERTEX_ARRAY);
    }

    // Draw range query box
    if (m_drawRangeQueryBox)
    {
        // Draw range box
        glPushAttrib(GL_POLYGON_BIT);
        glColor4ubv(m_RQ_box_color);
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

    // Draw range query result
    if (m_drawRangeQueryResult && !m_pointsInRange.empty())
    {
        // Draw points
        glEnableClientState(GL_VERTEX_ARRAY);

        glPointSize(m_RQ_result_size);
        glColor4ubv(m_RQ_result_color);
        glVertexPointer(3, GL_DOUBLE, sizeof(Point3d), &m_pointsInRange[0]);
        glDrawArrays(GL_POINTS, 0, (unsigned int)m_pointsInRange.size());

        glPointSize(m_PC_size);
        glDisableClientState(GL_VERTEX_ARRAY);
    }

    if (m_drawNearestQueryPoint)
    {
        glPointSize(m_NN_query_size);
        glColor4ubv(m_NN_query_color);

        glBegin(GL_POINTS);
            glVertex3dv(&m_nearestQueryPoint[0]);
        glEnd();
    }

    if (m_drawNearestResultPoint)
    {
        glPointSize(m_NN_result_size);
        glColor4ubv(m_NN_result_color);

        glBegin(GL_POINTS);
            glVertex3dv(&m_nearestResultPoint[0]);
        glEnd();
    }

    if (!m_smoothedPoints.empty() && m_drawSmoothedPoints)
    {
        glEnableClientState(GL_VERTEX_ARRAY);

        glPointSize(m_SM_size);

        if (m_drawSmoothedPointCloudWithColorArray)
        {
            glEnableClientState(GL_COLOR_ARRAY);
            glVertexPointer(3, GL_DOUBLE, sizeof(Point3d), &m_smoothedPoints[0]);
            glColorPointer(3, GL_UNSIGNED_BYTE, 3*sizeof(unsigned char), &m_pointColors[0]);
            glDrawArrays(GL_POINTS, 0, (unsigned int) m_smoothedPoints.size());

            glDisableClientState(GL_COLOR_ARRAY);
        }
        else
        {
            glColor4ubv(m_SM_color);
            glVertexPointer(3, GL_DOUBLE, sizeof(Point3d), &m_smoothedPoints[0]);
            glDrawArrays(GL_POINTS, 0, (unsigned int) m_smoothedPoints.size());
        }

        glDisableClientState(GL_VERTEX_ARRAY);
    }

    if (!m_thinnedPoints.empty() && m_drawThinnedPoints)
    {
        glEnableClientState(GL_VERTEX_ARRAY);

        glPointSize(m_TH_size);
        glColor4ubv(m_TH_color);
        glVertexPointer(3, GL_DOUBLE, sizeof(Point3d), &m_thinnedPoints[0]);
        glDrawArrays(GL_POINTS, 0, (unsigned int) m_thinnedPoints.size());

        glDisableClientState(GL_VERTEX_ARRAY);
    }

    if (m_drawTemporary)
    {
        glPointSize(m_TMPR_size);
        glColor4ubv(m_TMPS_color);

        glBegin(GL_POINTS);
            glVertex3dv(&m_tempPoint[0]);
        glEnd();

        if (!m_tempRadiusPoints.empty())
        {
            glEnableClientState(GL_VERTEX_ARRAY);

            glPointSize(m_TMPR_size);
            glColor4ubv(m_TMPR_color);
            glVertexPointer(3, GL_DOUBLE, sizeof(Point3d), &m_tempRadiusPoints[0]);
            glDrawArrays(GL_POINTS, 0, (unsigned int) m_tempRadiusPoints.size());

            glDisableClientState(GL_VERTEX_ARRAY);
        }
    }

    if (m_drawBFP)
    {
        glPushAttrib(GL_LINE_BIT);
        glLineWidth(2);
        glColor3f(0.0f, 1.0f, 0.0f);

        glBegin(GL_LINE_LOOP);
        for (std::size_t i = 0; i < 4; ++i)
        {
            glVertex3d(m_bfpCorners[i][0], m_bfpCorners[i][1], m_bfpCorners[i][2]);
        }
        glEnd();

        glPointSize(8);
        glBegin(GL_POINTS);
        for (std::size_t i = 0; i < 4; ++i)
        {
            glVertex3d(m_bfpCorners[i][0], m_bfpCorners[i][1], m_bfpCorners[i][2]);
        }
        glEnd();

		glPopAttrib();

        glColor4ub(129, 190, 255, 80);
        glBegin(GL_TRIANGLES);

        glVertex3d(m_bfpCorners[0][0], m_bfpCorners[0][1], m_bfpCorners[0][2]);
        glVertex3d(m_bfpCorners[1][0], m_bfpCorners[1][1], m_bfpCorners[1][2]);
        glVertex3d(m_bfpCorners[2][0], m_bfpCorners[2][1], m_bfpCorners[2][2]);

        glVertex3d(m_bfpCorners[2][0], m_bfpCorners[2][1], m_bfpCorners[2][2]);
        glVertex3d(m_bfpCorners[3][0], m_bfpCorners[3][1], m_bfpCorners[3][2]);
        glVertex3d(m_bfpCorners[0][0], m_bfpCorners[0][1], m_bfpCorners[0][2]);

        glEnd();
    }

	if (m_drawBFL)
	{
		glPushAttrib(GL_LINE_BIT);
		glLineWidth(4);
		glColor3f(79.0f / 255, 0.0f, 118.0f / 255);

		glBegin(GL_LINE_LOOP);
		for (std::size_t i = 0; i < 2; ++i)
		{
			glVertex3d(m_bflPoints[i][0], m_bflPoints[i][1], m_bflPoints[i][2]);
		}
		glEnd();
		glPopAttrib();
	}
    //draw coordinate frame
    drawCoordinateAxes();

    update();
}

bool GLwidget::eventFilter(QObject *obj, QEvent *event)
{
    if (event->type() == QEvent::KeyPress)
    {
        QKeyEvent* keyEvent = static_cast<QKeyEvent*>(event);

        if (keyEvent->key() == Qt::Key_H)
        {
            if (m_currentModeIsSmoothed)
            {
                m_drawSmoothedState = (m_drawSmoothedState + 1) % 6;

                auto minmax = std::minmax_element(m_distances.begin(), m_distances.end());
                double min = *minmax.first;
                double max = *minmax.second;

                switch (m_drawSmoothedState)
                {
                    case 0:
                        m_drawMainPointCloud = false;
                        m_drawSmoothedPointCloudWithColorArray = false;
                        m_drawSmoothedPoints = true;
                        break;
                    case 1:
                        m_drawMainPointCloud = true;
                        m_drawSmoothedPointCloudWithColorArray = false;
                        m_drawSmoothedPoints = true;
                        break;
                    case 2:
                        computeColorGrey(min, max);
                        m_drawMainPointCloud = false;
                        m_drawSmoothedPointCloudWithColorArray = true;
                        m_drawSmoothedPoints = true;
                        break;
                    case 3:
                        computeColorRainbow(min, max);
                        m_drawMainPointCloud = false;
                        m_drawSmoothedPointCloudWithColorArray = true;
                        m_drawSmoothedPoints = true;
                        break;
                    case 4:
                        computeColorHeat(min, max);
                        m_drawMainPointCloud = false;
                        m_drawSmoothedPointCloudWithColorArray = true;
                        m_drawSmoothedPoints = true;
                        break;
                    case 5:
                        m_drawMainPointCloud = true;
                        m_drawSmoothedPointCloudWithColorArray = false;
                        m_drawSmoothedPoints = false;
                        break;
                }
            }
            else
            {
                m_drawMainPointCloud = !m_drawMainPointCloud;
                m_drawSmoothedPointCloudWithColorArray = false;
            }

            this->update();
        }
    }

    return QObject::eventFilter(obj, event);
}

void GLwidget::reloadDrawSettings()
{
    loadDrawSettings();
    this->update();
}

void GLwidget::drawingMainPointCloudChange(bool value)
{
    m_drawMainPointCloud = value;
    this->update();
}

void GLwidget::drawingTemporaryChange(bool value)
{
    m_drawTemporary = value;
    this->update();
}

void GLwidget::drawingRangeQueryBoxChange(bool value)
{
    m_drawRangeQueryBox = value;

    if (!value)
        m_drawRangeQueryResult = value;

    this->update();
}

void GLwidget::drawingRangeQueryResultEnabled(bool value)
{
    m_drawRangeQueryResult = value;

    if (!value)
        m_pointsInRange.clear();

    this->update();
}

void GLwidget::rangeQueryCenterChanged(double x, double y, double z)
{
    m_rangeCenter[0] = x;
    m_rangeCenter[1] = y;
    m_rangeCenter[2] = z;

    updateRangeQueryBoxData();
    this->update();
}

void GLwidget::rangeQueryExtendChanged(double dx, double dy, double dz)
{
    m_rangeExtend[0] = dx;
    m_rangeExtend[1] = dy;
    m_rangeExtend[2] = dz;

    updateRangeQueryBoxData();
    this->update();
}

void GLwidget::rangeQueryResultChanged(std::vector<Point3d> points)
{
    m_pointsInRange = std::move(points);
    this->update();
}

void GLwidget::drawingNearestNeighborQueryPointChanged(bool value)
{
    m_drawNearestQueryPoint = value;

    if (!value)
        m_drawNearestResultPoint = value;

    this->update();
}

void GLwidget::drawingNearestNeighborResultPointChanged(bool value)
{
    m_drawNearestResultPoint = value;
    this->update();
}

void GLwidget::nearestNeighborQueryPointChanged(const Point3d& queryPoint)
{
    m_nearestQueryPoint = queryPoint;
    this->update();
}

void GLwidget::nearestNeighborResultPointChanged(const Point3d& resultPoint)
{
    m_nearestResultPoint = resultPoint;
    this->update();
}

void GLwidget::drawingSmoothedPointsChanged(bool value)
{
    m_drawSmoothedPoints = value;
    m_currentModeIsSmoothed = value;
    this->update();
}

void GLwidget::drawingThinnedPointsChanged(bool value)
{
    m_drawThinnedPoints = value;
    this->update();
}

void GLwidget::drawingBestFitPlaneChanged(bool value)
{
    m_drawBFP = value;
	m_drawBFL = value;
    this->update();
}

void GLwidget::loadDrawSettings()
{
    QFileInfo check_file(QString("drawSettings.ini"));

    if (!(check_file.exists() && check_file.isFile()))
        std::cerr << "Could not find drawSettings.ini file in executable directory."
                     "Use default values!\n";

    QSettings set(QString("drawSettings.ini"), QSettings::IniFormat);

    m_SM_color[0] = set.value("SM_color_R", 255).value<unsigned char>();
    m_SM_color[1] = set.value("SM_color_G", 0).value<unsigned char>();
    m_SM_color[2] = set.value("SM_color_B", 0).value<unsigned char>();
    m_SM_color[3] = set.value("SM_color_A", 255).value<unsigned char>();
    m_SM_size = set.value("SM_size", 2).toInt();

    m_TH_color[0] = set.value("TH_color_R", 0).value<unsigned char>();
    m_TH_color[1] = set.value("TH_color_G", 255).value<unsigned char>();
    m_TH_color[2] = set.value("TH_color_B", 0).value<unsigned char>();
    m_TH_color[3] = set.value("TH_color_A", 255).value<unsigned char>();
    m_TH_size = set.value("TH_size", 2).toInt();

    m_TMPS_color[0] = set.value("TMPS_color_R", 255).value<unsigned char>();
    m_TMPS_color[1] = set.value("TMPS_color_G", 255).value<unsigned char>();
    m_TMPS_color[2] = set.value("TMPS_color_B", 255).value<unsigned char>();
    m_TMPS_color[3] = set.value("TMPS_color_A", 255).value<unsigned char>();
    m_TMPS_size = set.value("TMPS_size", 4).toInt();

    m_TMPR_color[0] = set.value("TMPR_color_R", 0).value<unsigned char>();
    m_TMPR_color[1] = set.value("TMPR_color_G", 0).value<unsigned char>();
    m_TMPR_color[2] = set.value("TMPR_color_B", 0).value<unsigned char>();
    m_TMPR_color[3] = set.value("TMPR_color_A", 255).value<unsigned char>();
    m_TMPR_size = set.value("TMPR_size", 4).toInt();

    m_PC_color[0] = set.value("PC_color_R", 255).value<unsigned char>();
    m_PC_color[1] = set.value("PC_color_G", 141).value<unsigned char>();
    m_PC_color[2] = set.value("PC_color_B", 42).value<unsigned char>();
    m_PC_color[3] = set.value("PC_color_A", 255).value<unsigned char>();
    m_PC_size = set.value("PC_size", 2).toInt();

    m_RQ_box_color[0] = set.value("RQ_box_color_R", 188).value<unsigned char>();
    m_RQ_box_color[1] = set.value("RQ_box_color_G", 217).value<unsigned char>();
    m_RQ_box_color[2] = set.value("RQ_box_color_B", 5).value<unsigned char>();
    m_RQ_box_color[3] = set.value("RQ_box_color_A", 255).value<unsigned char>();
    m_RQ_box_width = set.value("RQ_box_width", 1.0).toDouble();

    m_RQ_result_color[0] = set.value("RQ_result_color_R", 217).value<unsigned char>();
    m_RQ_result_color[1] = set.value("RQ_result_color_G", 22).value<unsigned char>();
    m_RQ_result_color[2] = set.value("RQ_result_color_B", 25).value<unsigned char>();
    m_RQ_result_color[3] = set.value("RQ_result_color_A", 255).value<unsigned char>();
    m_RQ_result_size = set.value("RQ_result_size", 2).toInt();

    m_NN_query_color[0] = set.value("NN_query_color_R", 255).value<unsigned char>();
    m_NN_query_color[1] = set.value("NN_query_color_G", 0).value<unsigned char>();
    m_NN_query_color[2] = set.value("NN_query_color_B", 255).value<unsigned char>();
    m_NN_query_color[3] = set.value("NN_query_color_A", 255).value<unsigned char>();
    m_NN_query_size = set.value("NN_query_size", 15).toInt();

    m_NN_result_color[0] = set.value("NN_result_color_R", 0).value<unsigned char>();
    m_NN_result_color[1] = set.value("NN_result_color_G", 0).value<unsigned char>();
    m_NN_result_color[2] = set.value("NN_result_color_B", 255).value<unsigned char>();
    m_NN_result_color[3] = set.value("NN_result_color_A", 255).value<unsigned char>();
    m_NN_result_size = set.value("NN_result_size", 15).toInt();

    m_XAXIS_color[0] = set.value("XAXIS_color_R", 255).value<unsigned char>();
    m_XAXIS_color[1] = set.value("XAXIS_color_G", 0).value<unsigned char>();
    m_XAXIS_color[2] = set.value("XAXIS_color_B", 0).value<unsigned char>();
    m_XAXIS_color[3] = set.value("XAXIS_color_A", 255).value<unsigned char>();

    m_YAXIS_color[0] = set.value("YAXIS_color_R", 0).value<unsigned char>();
    m_YAXIS_color[1] = set.value("YAXIS_color_G", 255).value<unsigned char>();
    m_YAXIS_color[2] = set.value("YAXIS_color_B", 0).value<unsigned char>();
    m_YAXIS_color[3] = set.value("YAXIS_color_A", 255).value<unsigned char>();

    m_ZAXIS_color[0] = set.value("ZAXIS_color_R", 0).value<unsigned char>();
    m_ZAXIS_color[1] = set.value("ZAXIS_color_G", 0).value<unsigned char>();
    m_ZAXIS_color[2] = set.value("ZAXIS_color_B", 255).value<unsigned char>();
    m_ZAXIS_color[3] = set.value("ZAXIS_color_A", 255).value<unsigned char>();
    m_AXIS_width = set.value("AXIS_width", 1.0).toDouble();

    m_CENTERSPHERE_color[0] = set.value("CENTERSPHERE_color_R", 255).value<unsigned char>();
    m_CENTERSPHERE_color[1] = set.value("CENTERSPHERE_color_G", 255).value<unsigned char>();
    m_CENTERSPHERE_color[2] = set.value("CENTERSPHERE_color_B", 0).value<unsigned char>();
    m_CENTERSPHERE_color[3] = set.value("CENTERSPHERE_color_A", 255).value<unsigned char>();

    m_XCIRCLE_color[0] = set.value("XCIRCLE_color_R", 255).value<unsigned char>();
    m_XCIRCLE_color[1] = set.value("XCIRCLE_color_G", 0).value<unsigned char>();
    m_XCIRCLE_color[2] = set.value("XCIRCLE_color_B", 0).value<unsigned char>();
    m_XCIRCLE_color[3] = set.value("XCIRCLE_color_A", 255).value<unsigned char>();

    m_YCIRCLE_color[0] = set.value("YCIRCLE_color_R", 0).value<unsigned char>();
    m_YCIRCLE_color[1] = set.value("YCIRCLE_color_G", 255).value<unsigned char>();
    m_YCIRCLE_color[2] = set.value("YCIRCLE_color_B", 0).value<unsigned char>();
    m_YCIRCLE_color[3] = set.value("YCIRCLE_color_A", 255).value<unsigned char>();

    m_ZCIRCLE_color[0] = set.value("ZCIRCLE_color_R", 0).value<unsigned char>();
    m_ZCIRCLE_color[1] = set.value("ZCIRCLE_color_G", 0).value<unsigned char>();
    m_ZCIRCLE_color[2] = set.value("ZCIRCLE_color_B", 255).value<unsigned char>();
    m_ZCIRCLE_color[3] = set.value("ZCIRCLE_color_A", 255).value<unsigned char>();

    m_BB_color[0] = set.value("BB_color_R", 255).value<unsigned char>();
    m_BB_color[1] = set.value("BB_color_G", 255).value<unsigned char>();
    m_BB_color[2] = set.value("BB_color_B", 255).value<unsigned char>();
    m_BB_color[3] = set.value("BB_color_A", 255).value<unsigned char>();
    m_BB_width = set.value("BB_width", 1.0).toDouble();

    m_BG_top_color[0] = set.value("BG_top_color_R", 38).value<unsigned char>();
    m_BG_top_color[1] = set.value("BG_top_color_G", 38).value<unsigned char>();
    m_BG_top_color[2] = set.value("BG_top_color_B", 38).value<unsigned char>();
    m_BG_top_color[3] = set.value("BG_top_color_A", 255).value<unsigned char>();

    m_BG_bottom_color[0] = set.value("BG_bottom_color_R", 100).value<unsigned char>();
    m_BG_bottom_color[1] = set.value("BG_bottom_color_G", 100).value<unsigned char>();
    m_BG_bottom_color[2] = set.value("BG_bottom_color_B", 100).value<unsigned char>();
    m_BG_bottom_color[3] = set.value("BG_bottom_color_A", 255).value<unsigned char>();
}

void GLwidget::mousePressEvent(QMouseEvent * e)
{
    if (e->buttons() == Qt::LeftButton)
        m_mouseLastPos = e->pos();
}

void GLwidget::mouseMoveEvent(QMouseEvent * e)
{
    if (e->buttons() != Qt::LeftButton){ return; }

    //No action, if mouse position outside the window
    if ((*e).x() < 0 || (*e).x() >= width()){ return; }
    if ((*e).y() < 0 || (*e).y() >= height()){ return; }

    const QPoint& mouseCurrentPos(e->pos());

    if (m_mouseLastPos == mouseCurrentPos) return;

    makeCurrent();
    m_camera.rotate(m_mouseLastPos.x(), m_mouseLastPos.y(),
                    mouseCurrentPos.x(), mouseCurrentPos.y());

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
        const Point3d& pt = m_points[i];
        if (pt[0] < minPoint[0]) minPoint[0] = pt[0];
        else if (pt[0] > maxPoint[0]) maxPoint[0] = pt[0];

        if (pt[1] < minPoint[1]) minPoint[1] = pt[1];
        else if (pt[1] > maxPoint[1]) maxPoint[1] = pt[1];

        if (pt[2] < minPoint[2]) minPoint[2] = pt[2];
        else if (pt[2] > maxPoint[2]) maxPoint[2] = pt[2];
    }

    m_sceneCenter[0] = (maxPoint[0] + minPoint[0]) / 2;
    m_sceneCenter[1] = (maxPoint[1] + minPoint[1]) / 2;
    m_sceneCenter[2] = (maxPoint[2] + minPoint[2]) / 2;

    m_sceneRadius = distance3d(m_sceneCenter, maxPoint);

    m_bbmin=minPoint;
    m_bbmax=maxPoint;

    makeCurrent();
    //set rotation center and scene radius to initially setup the camera
    m_camera.initializeCamera(m_sceneCenter, m_sceneRadius);
    //updateProjection to make the new point cloud fit to the screen
    m_camera.updateProjection();
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
    glColor4ubv(m_XAXIS_color);
    glVertex3dv(&m_sceneCenter[0]);
    glVertex3d(m_sceneCenter[0] + m_sceneRadius, m_sceneCenter[1], m_sceneCenter[2]);
    //draw line for Y-Axis
    glColor4ubv(m_YAXIS_color);
    glVertex3dv(&m_sceneCenter[0]);
    glVertex3d(m_sceneCenter[0], m_sceneCenter[1] + m_sceneRadius, m_sceneCenter[2]);
    //draw line for Z-Axis
    glColor4ubv(m_ZAXIS_color);
    glVertex3dv(&m_sceneCenter[0]);
    glVertex3d(m_sceneCenter[0], m_sceneCenter[1], m_sceneCenter[2] + m_sceneRadius);
    glEnd();

    //draw center point as a sphere
    glPushMatrix();
    glColor4ubv(m_CENTERSPHERE_color);
    glTranslated(m_sceneCenter[0], m_sceneCenter[1], m_sceneCenter[2]);
    GLUquadric* quad = gluNewQuadric();
    gluSphere(quad, m_sceneRadius / 20, 30, 30);
    gluDeleteQuadric(quad);
    glPopMatrix();

    glPushMatrix();
    glColor4ubv(m_ZCIRCLE_color);
    glTranslated(m_sceneCenter[0], m_sceneCenter[1], m_sceneCenter[2]);
    glScaled(m_sceneRadius, m_sceneRadius, m_sceneRadius);
    drawCircle();
    //draw another circle 90 degree rotated
    glRotated(90, 1, 0, 0);
    glColor4ubv(m_YCIRCLE_color);
    drawCircle();
    glRotated(90, 0, 1, 0);
    glColor4ubv(m_XCIRCLE_color);
    drawCircle();
    glPopMatrix();

    //draw bounding box
    glPushMatrix();
    glPushAttrib(GL_POLYGON_BIT);
    glColor4ubv(m_BB_color);
    Point3d S=m_bbmax-m_bbmin;
    glTranslated(m_bbmin[0], m_bbmin[1], m_bbmin[2]);
    glScaled(S[0], S[1], S[2]);
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
    glColor4ubv(m_BG_bottom_color); //color bottom
    glVertex2f(0.0f, 0.0f);  glVertex2f(winWidth, 0.0f);
    glColor4ubv(m_BG_top_color);  //color top
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

    m_fbl[0] = m_rangeCenter[0] - (m_rangeExtend[0] / 2.0);
    m_fbl[1] = m_rangeCenter[1] - (m_rangeExtend[1] / 2.0);
    m_fbl[2] = m_rangeCenter[2] - (m_rangeExtend[2] / 2.0);

    m_fbr[0] = m_rangeCenter[0] + (m_rangeExtend[0] / 2.0);
    m_fbr[1] = m_rangeCenter[1] - (m_rangeExtend[1] / 2.0);
    m_fbr[2] = m_rangeCenter[2] - (m_rangeExtend[2] / 2.0);

    m_ftl[0] = m_rangeCenter[0] - (m_rangeExtend[0] / 2.0);
    m_ftl[1] = m_rangeCenter[1] - (m_rangeExtend[1] / 2.0);
    m_ftl[2] = m_rangeCenter[2] + (m_rangeExtend[2] / 2.0);

    m_ftr[0] = m_rangeCenter[0] + (m_rangeExtend[0] / 2.0);
    m_ftr[1] = m_rangeCenter[1] - (m_rangeExtend[1] / 2.0);
    m_ftr[2] = m_rangeCenter[2] + (m_rangeExtend[2] / 2.0);

    // Back

    m_bbl[0] = m_rangeCenter[0] - (m_rangeExtend[0] / 2.0);
    m_bbl[1] = m_rangeCenter[1] + (m_rangeExtend[1] / 2.0);
    m_bbl[2] = m_rangeCenter[2] - (m_rangeExtend[2] / 2.0);

    m_bbr[0] = m_rangeCenter[0] + (m_rangeExtend[0] / 2.0);
    m_bbr[1] = m_rangeCenter[1] + (m_rangeExtend[1] / 2.0);
    m_bbr[2] = m_rangeCenter[2] - (m_rangeExtend[2] / 2.0);

    m_btl[0] = m_rangeCenter[0] - (m_rangeExtend[0] / 2.0);
    m_btl[1] = m_rangeCenter[1] + (m_rangeExtend[1] / 2.0);
    m_btl[2] = m_rangeCenter[2] + (m_rangeExtend[2] / 2.0);

    m_btr[0] = m_rangeCenter[0] + (m_rangeExtend[0] / 2.0);
    m_btr[1] = m_rangeCenter[1] + (m_rangeExtend[1] / 2.0);
    m_btr[2] = m_rangeCenter[2] + (m_rangeExtend[2] / 2.0);
}

void GLwidget::computeColorGrey(double min, double max)
{
#pragma omp parallel for
    for (int i = 0; i < m_distances.size(); ++i)
    {
        unsigned char factor = ((m_distances[i] - min) / (max - min)) * 255;
        m_pointColors[i*3] = factor;
        m_pointColors[i*3 + 1] = factor;
        m_pointColors[i*3 + 2] = factor;
    }
}

void GLwidget::computeColorRainbow(double min, double max)
{
    const std::vector<std::array<unsigned char, 3>> lookup {{0,0,255},
                                                            {0,255,0},
                                                            {255,255,0},
                                                            {255,127,0},
                                                            {255,0,0}};


    const double spacing = 1.0 / (lookup.size() - 1);

#pragma omp parallel for
    for (int i = 0; i < m_distances.size(); ++i)
    {
        double factor = ((m_distances[i] - min) / (max - min));

        std::size_t startIdx = std::floor((lookup.size() - 1) * factor);
        std::size_t endIdx = startIdx + 1;

        double startFactor = startIdx * spacing;
        double endFactor = endIdx * spacing;

        double factorBetween = (factor - startFactor) / (endFactor - startFactor);

        m_pointColors[i*3] = (1.0 - factorBetween)*lookup[startIdx][0] + factorBetween*lookup[endIdx][0];
        m_pointColors[i*3 + 1] = (1.0 - factorBetween)*lookup[startIdx][1] + factorBetween*lookup[endIdx][1];
        m_pointColors[i*3 + 2] = (1.0 - factorBetween)*lookup[startIdx][2] + factorBetween*lookup[endIdx][2];
    }
}

void GLwidget::computeColorHeat(double min, double max)
{
	const std::vector<std::array<unsigned char, 3>> lookup {{0,0,0},
															{255,0,0},
															{255,255,0},
															{255,255,255}} ;

    const double spacing = 1.0 / (lookup.size() - 1);

#pragma omp parallel for
    for (int i = 0; i < m_distances.size(); ++i)
    {
        double factor = ((m_distances[i] - min) / (max - min));

        std::size_t startIdx = std::floor((lookup.size() - 1) * factor);
        std::size_t endIdx = startIdx + 1;

        double startFactor = startIdx * spacing;
        double endFactor = endIdx * spacing;

        double factorBetween = (factor - startFactor) / (endFactor - startFactor);

        m_pointColors[i*3] = (1.0 - factorBetween)*lookup[startIdx][0] + factorBetween*lookup[endIdx][0];
        m_pointColors[i*3 + 1] = (1.0 - factorBetween)*lookup[startIdx][1] + factorBetween*lookup[endIdx][1];
        m_pointColors[i*3 + 2] = (1.0 - factorBetween)*lookup[startIdx][2] + factorBetween*lookup[endIdx][2];
    }
}

void GLwidget::computeColorDiverge(double min, double max)
{
    const std::vector<std::array<unsigned char, 3>> lookup {{142, 1, 82},
                                                            {247, 247, 247},
                                                            {39, 100, 25}};

#pragma omp parallel for
    for (int i = 0; i < m_distances.size(); ++i)
    {
        if (m_distances[i] <= 0)
        {
            const double factor = ((m_distances[i] - min) / (-min));

            m_pointColors[i*3] = (1.0 - factor)*lookup[0][0] + factor*lookup[1][0];
            m_pointColors[i*3 + 1] = (1.0 - factor)*lookup[0][1] + factor*lookup[1][1];
            m_pointColors[i*3 + 2] = (1.0 - factor)*lookup[0][2] + factor*lookup[1][2];
        }
        else
        {
            const double factor = ((m_distances[i]) / (max));

            m_pointColors[i*3] = (1.0 - factor)*lookup[1][0] + factor*lookup[2][0];
            m_pointColors[i*3 + 1] = (1.0 - factor)*lookup[1][1] + factor*lookup[2][1];
            m_pointColors[i*3 + 2] = (1.0 - factor)*lookup[1][2] + factor*lookup[2][2];
        }
    }
}
