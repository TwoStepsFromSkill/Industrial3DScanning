#ifndef MY_GLwidget_H
#define MY_GLwidget_H

#include <QtWidgets/QOpenGLWidget>
#include <string>     //we want to process text + filenames
#include <iostream>   //for making output text to the console with "cout"
#include <vector>     //for std::vector class functions
#include <stdio.h>

#include "Point3d.h"
#include "GLcamera.h"

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
    void setPoints    (std::vector<Point3d>& points)   { m_points = points; updateScene(); }

    //access to data
    std::vector<Point3d>& points() { return m_points; }

    void setSmoothedPoints(const std::vector<Point3d>& smoothedPoints) { m_smoothedPoints = smoothedPoints; }

    //return camera
    GLcamera& camera(){return m_camera;}

public slots:
    // Visualization
    void reloadDrawSettings();

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

private:
    void loadDrawSettings();
    void writeSettings();

    void  mousePressEvent(QMouseEvent * e);  ///<
    void  mouseMoveEvent (QMouseEvent * e);  ///<
    void  wheelEvent     (QWheelEvent * e);  ///<

    void drawBox();             ///< draws a unit box
    void drawCircle();          ///< draws a unit circle

    void drawCoordinateAxes();  ///< draws the coordinate system
    void drawBackground();      ///< draws the scene background

    void updateRangeQueryBoxData();

    std::vector<Point3d> m_points;    //point data
    std::vector<Point3d> m_pointsInRange;
    std::vector<Point3d> m_smoothedPoints;
    std::vector<Point3d> m_thinnedPoints;

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

    // Nearest Neighbor
    Point3d m_nearestQueryPoint;
    Point3d m_nearestResultPoint;
    bool m_drawNearestQueryPoint;
    bool m_drawNearestResultPoint;

    // Smoothing
    bool m_drawSmoothedPoints;

    // Colors and draw settings
    unsigned char m_PC_color[4];
    int m_PC_size;

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
