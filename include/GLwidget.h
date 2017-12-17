#ifndef MY_GLwidget_H
#define MY_GLwidget_H

#include <QtWidgets/QOpenGLWidget>
#include <string>     //we want to process text + filenames
#include <iostream>   //for making output text to the console with "cout"
#include <vector>     //for std::vector class functions
#include <stdio.h>

#include "Point3d.h"
#include "GLcamera.h"

enum class ColorScale
{
    GREY,
    RAINBOW,
    HEAT,
    DIVERGE
};

class GLwidget : public QOpenGLWidget
{
    Q_OBJECT
public:
    GLwidget(QWidget* parent = nullptr);

    //overloaded QT functions
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();

    //updates size of the scene
    void updateScene();

    //points to draw
    void setPoints    (std::vector<Point3d>& points)
    {
        m_points = points;
        m_pointColors.resize(points.size() * 3);
        updateScene();
    }
    void setPointDistances (std::vector<double>& distances) { m_distances = distances; }

    //access to data
    std::vector<Point3d>& points() { return m_points; }

    void setSmoothedPoints(const std::vector<Point3d>& smoothedPoints) { m_smoothedPoints = smoothedPoints; }
    void setThinnedPoints(const std::vector<Point3d>& thinnedPoints) { m_thinnedPoints = thinnedPoints; }

    void setBFPCorners(const std::vector<Point3d>& corners) { m_bfpCorners = corners; }

	void setBFLPoints(const std::vector<Point3d>& points) { m_bflPoints = points; }

	void setBFSPoints(const std::vector<Point3d>& points) { m_spherePoints = points; }

    void setTempPoint(const Point3d& p) { m_tempPoint = p; }
    void setTempRadiusPoints(const std::vector<Point3d>& pts) { m_tempRadiusPoints = pts; }

    void drawingMainCloudPointWithColorArray(bool val)
    {
        m_drawMainPointCloudWithColorArray = val;
        updateColors();
        this->update();
    }

    void updateColors()
    {
        auto minmax = std::minmax_element(m_distances.begin(), m_distances.end());
        double min = *minmax.first;
        double max = *minmax.second;

        switch (m_colorScale)
        {
            case ColorScale::GREY:
                computeColorDiverge(min, max);
                break;
            case ColorScale::RAINBOW:
                computeColorRainbow(min, max);
                break;
            case ColorScale::HEAT:
                computeColorHeat(min, max);
                break;
            case ColorScale::DIVERGE:
                computeColorDiverge(min, max);
                break;
        }
    }

    void colorScaleToGrey()
    {
        m_colorScale = ColorScale::GREY;
        updateColors();
    }
    void colorScaleToRainbow()
    {
        m_colorScale = ColorScale::RAINBOW;
        updateColors();
    }
    void colorScaleToHeat()
    {
        m_colorScale = ColorScale::HEAT;
        updateColors();
    }
    void colorScaleToDiverge()
   {
        m_colorScale = ColorScale::DIVERGE;
        updateColors();
    }

    void drawingSmoothedPointWithColorArray(bool val)
    { m_drawSmoothedPointCloudWithColorArray = val; this->update(); }

    //return camera
    GLcamera& camera(){return m_camera;}

protected:
    bool eventFilter(QObject *obj, QEvent *event);

public slots:
    // Visualization
    void reloadDrawSettings();

    void drawingMainPointCloudChange(bool value);
    void drawingTemporaryChange(bool value);

    // Range query
    void drawingRangeQueryBoxChange(bool value);
    void drawingRangeQueryResultEnabled(bool value);
    void rangeQueryCenterChanged(double x, double y, double z);
    void rangeQueryExtendChanged(double dx, double dy, double dz);
    void rangeQueryResultChanged(std::vector<Point3d> points);

    // Nearest neighbor
    void drawingNearestNeighborQueryPointChanged(bool value);
    void drawingNearestNeighborResultPointChanged(bool value);
    void nearestNeighborQueryPointChanged(const Point3d& queryPoint);
    void nearestNeighborResultPointChanged(const Point3d& resultPoint);

    // Smoothed points
    void drawingSmoothedPointsChanged(bool value);

    // Thinned points
    void drawingThinnedPointsChanged(bool value);

    void drawingBestFitLineChanged(bool value);
    void drawingBestFitPlaneChanged(bool value);
    void drawingBestFitSphereChanged(bool value);

private:
    void loadDrawSettings();
    void writeSettings();

    void mousePressEvent(QMouseEvent * e);
    void mouseMoveEvent (QMouseEvent * e);
    void wheelEvent     (QWheelEvent * e);

    void drawBox();             ///< draws a unit box
    void drawCircle();          ///< draws a unit circle

    void drawCoordinateAxes();  ///< draws the coordinate system
    void drawBackground();      ///< draws the scene background

    void updateRangeQueryBoxData();

    void computeColorGrey(double min, double max);
    void computeColorRainbow(double min, double max);
    void computeColorHeat(double min, double max);
    void computeColorDiverge(double min, double max);

    ColorScale m_colorScale;

    std::vector<Point3d> m_points;    //point data
    std::vector<unsigned char> m_pointColors;
    std::vector<double> m_distances;
    std::vector<Point3d> m_pointsInRange;
    std::vector<Point3d> m_smoothedPoints;
    std::vector<Point3d> m_thinnedPoints;

    std::vector<Point3d> m_bfpCorners;
    bool m_drawBFP;

	std::vector<Point3d> m_bflPoints;
	bool m_drawBFL;

    std::vector<Point3d> m_spherePoints;
    bool m_drawBFS;

    Point3d m_tempPoint;
    std::vector<Point3d> m_tempRadiusPoints;

    QPoint               m_mouseLastPos;  //last mouse position clicked

    GLcamera  m_camera;         //virtual camera
    Point3d   m_bbmin,m_bbmax;  //bounding box coordinates
    Point3d   m_sceneCenter;    //center of the scene
    double    m_sceneRadius;    //radius of the scene

    // Range query
    Point3d m_rangeCenter;
    Point3d m_rangeExtend;
    bool m_drawRangeQueryBox;
    bool m_drawRangeQueryResult;

    double m_fbl[3];
    double m_fbr[3];
    double m_ftl[3];
    double m_ftr[3];

    double m_bbl[3];
    double m_bbr[3];
    double m_btl[3];
    double m_btr[3];

    // Point cloud
    bool m_drawMainPointCloud;
    bool m_drawMainPointCloudWithColorArray;
    bool m_drawSmoothedPointCloudWithColorArray;
    std::size_t m_drawSmoothedState;

    bool m_drawTemporary;

    // Nearest Neighbor
    Point3d m_nearestQueryPoint;
    Point3d m_nearestResultPoint;
    bool m_drawNearestQueryPoint;
    bool m_drawNearestResultPoint;

    // Smoothing
    bool m_currentModeIsSmoothed;
    bool m_drawSmoothedPoints;
    unsigned char m_SM_color[4];
    int m_SM_size;

    // Thinning
    bool m_drawThinnedPoints;
    unsigned char m_TH_color[4];
    int m_TH_size;

    // Colors and draw settings
    unsigned char m_PC_color[4];
    int m_PC_size;

    // Temp
    unsigned char m_TMPS_color[4];
    int m_TMPS_size;

    unsigned char m_TMPR_color[4];
    int m_TMPR_size;

    unsigned char m_RQ_box_color[4];
    double m_RQ_box_width;

    unsigned char m_RQ_result_color[4];
    int m_RQ_result_size;

    unsigned char m_NN_query_color[4];
    int m_NN_query_size;

    unsigned char m_NN_result_color[4];
    int m_NN_result_size;

    unsigned char m_XAXIS_color[4];
    unsigned char m_YAXIS_color[4];
    unsigned char m_ZAXIS_color[4];
    double m_AXIS_width;

    unsigned char m_CENTERSPHERE_color[4];

    unsigned char m_XCIRCLE_color[4];
    unsigned char m_YCIRCLE_color[4];
    unsigned char m_ZCIRCLE_color[4];

    unsigned char m_BB_color[4];
    double m_BB_width;

    unsigned char m_BG_top_color[4];
    unsigned char m_BG_bottom_color[4];
};

#endif
