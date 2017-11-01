#ifndef MY_MAINWINDOW_H
#define MY_MAINWINDOW_H

#include <QtWidgets/QMainWindow>
#include "GLwidget.h"
#include "Point3d.h"
#include "KDTree.h"
#include "RangeQueryDialog.h"

#include <vector>

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow();

private:
    QMenu*    m_fileMenu;
    QMenu*    m_viewMenu;
    QMenu*    m_rangeMenu;
    RangeQueryDialog* m_rangeDialog;
    GLwidget* m_glWidget;

    std::vector<Point3d> m_points;

private:
    void loadFileXYZ(const char* filename, std::vector<Point3d>& points);

private slots:
    void openFile();
	std::vector<Point3d> getPointsInRange(Node* KDTree, Point3d* min, Point3d* max);
	void searchKDTree(Node * curNode, unsigned int depth, Point3d * rMin, Point3d * rMax, std::vector<Point3d>* queryVector);
	void changeProjection();
    void setRangeQuery();
};

#endif
