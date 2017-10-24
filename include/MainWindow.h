#ifndef MY_MAINWINDOW_H
#define MY_MAINWINDOW_H

#include <QtWidgets/QMainWindow>
#include "GLwidget.h"
#include "Point3d.h"
#include "KDTree.h"

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow();

private:
    QMenu*    m_fileMenu;
    QMenu*    m_viewMenu;
    GLwidget* m_glWidget;

private:
    void loadFileXYZ(const char* filename, std::vector<Point3d>& points);

private slots:
    void openFile();
    void changeProjection();
};

#endif
