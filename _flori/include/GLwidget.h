#ifndef MY_GLwidget_H
#define MY_GLwidget_H

#include <QtWidgets\QOpenGLWidget>
#include <string>     //we want to process text + filenames
#include <iostream>   //for making output text to the console with "cout"
#include <vector>     //for std::vector class functions
#include <stdio.h>

#include "Point3d.h"
#include "GLcamera.h"

class GLwidget : public QOpenGLWidget
{
  public:
    //overloaded QT functions
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();

    //updates size of the scene
    void updateScene();
    
    //points to draw
    std::vector<Point3d>& points(){return m_points;}

    //return camera 
    GLcamera& camera(){return m_camera;}

  private:
    void  mousePressEvent(QMouseEvent * e);  ///< 
    void  mouseMoveEvent (QMouseEvent * e);  ///< 
    void  wheelEvent     (QWheelEvent * e);  ///<

    void drawBox();

    std::vector<Point3d> m_points;        //point data
    QPoint               m_mouseLastPos;  //last mouse position clicked
    
    GLcamera  m_camera;         //virtual camera
    Point3d   m_bbmin,m_bbmax;  //bounding box coordinates
    Point3d   m_sceneCenter;    //center of the scene
    double    m_sceneRadius;    //radius of the scene
};

#endif