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
    std::vector<Point3d>& points() { return m_points; } //return reference not the copy!

    //return camera
    GLcamera& camera(){return m_camera;}

    void setPointsInRange(const std::vector<Point3d>& points);

public slots:
    void drawingRangeQueryBoxChange(bool value);
    void drawingRangeQueryResultEnabled(bool value);
    void rangeQueryCenterChanged(double x, double y, double z);
    void rangeQueryExtendChanged(double dx, double dy, double dz);

private:
    void  mousePressEvent(QMouseEvent * e);  ///<
    void  mouseMoveEvent (QMouseEvent * e);  ///<
    void  wheelEvent     (QWheelEvent * e);  ///<

    void drawBox();             ///< draws a unit box
    void drawCircle();          ///< draws a unit circle

    void drawCoordinateAxes();  ///< draws the coordinate system
    void drawBackground();      ///< draws the scene background

    std::vector<Point3d> m_points;    //point data

    QPoint               m_mouseLastPos;  //last mouse position clicked

    GLcamera  m_camera;         //virtual camera
    Point3d   m_bbmin,m_bbmax;  //bounding box coordinates
    Point3d   m_sceneCenter;    //center of the scene
    double    m_sceneRadius;    //radius of the scene

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

    std::vector<Point3d> m_pointsInRange;

    void updateRangeQueryBoxData();
};

#endif
