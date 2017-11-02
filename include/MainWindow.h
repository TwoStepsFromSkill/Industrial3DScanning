#ifndef MY_MAINWINDOW_H
#define MY_MAINWINDOW_H

#include <QtWidgets/QMainWindow>
#include "GLwidget.h"
#include "Point3d.h"
#include "KDTree.h"

#include <vector>

class QHBoxLayout;
class QTabWidget;
class RangeQueryWidget;
class NearestNeighboorWidget;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow();

private:
    QMenu*    m_fileMenu;
    QMenu*    m_viewMenu;

    QHBoxLayout* m_mainLayout;

    GLwidget* m_glWidget;
    QTabWidget* m_tabWidget;

    RangeQueryWidget* m_rangeWidget;
    bool m_liveUpdateRangeQuery;

    NearestNeighboorWidget* m_nearestWidget;
    bool m_liveUpdateNearestNeighboor;

    std::vector<Point3d> m_points;
    Node* m_kdTree;

signals:
    void drawingRangeQueryResultChanged(bool);
    void rangeQueryCenterChange(double, double, double);
    void rangeQueryExtendChange(double, double, double);

private:
    void loadFileXYZ(const char* filename, std::vector<Point3d>& points);
    void updateSidebarWidgetData();

private slots:
    void openFile();
    void changeProjection();

    void tabSwitched(int);

    void rangeQueryCenterChanged(double, double, double);
    void rangeQueryExtendChanged(double, double, double);

    void rangeQueryLiveUpdateChange(bool value);

    void applyRangeQueryPressed();
    void hideRangeQueryPressed();

    void computeAndVisualizeRangeQuery();
};

#endif
