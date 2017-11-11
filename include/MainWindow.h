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
class NearestNeighborWidget;
class SmoothingWidget;
class ThinningWidget;

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

    NearestNeighborWidget* m_nearestWidget;
    bool m_liveUpdateNearestNeighbor;

    SmoothingWidget* m_smoothingWidget;

    ThinningWidget* m_thinningWidget;

    std::vector<Point3d> m_points;
    Node* m_kdTree;

signals:
    void reloadDrawSettings();

    void drawingRangeQueryResultChanged(bool);
    void rangeQueryCenterChange(double, double, double);
    void rangeQueryExtendChange(double, double, double);
    void rangeQueryResultChange(std::vector<Point3d>);

    void drawingNearestNeighborResultChanged(bool);
    void nearestNeighborQueryPointChange(const Point3d&);
    void nearestNeighborResultPointChange(const Point3d&);

private:
    void loadFileXYZ(const char* filename, std::vector<Point3d>& points);
    void updateSidebarWidgetData();

    void computeAndVisualizeRangeQuery();
	void testRangeRadius();
	void computeAndVisualizeNearestNeighbor();
	
private slots:
    void openFile();
    void changeProjection();

    void tabSwitched(int);

    // Range query
    void rangeQueryCenterChanged(double, double, double);
    void rangeQueryExtendChanged(double, double, double);

    void rangeQueryLiveUpdateChange(bool value);

    void applyRangeQueryPressed();
    void hideRangeQueryPressed();
	
    // Nearest neighbor
    void nearestNeighborQueryPointChanged(double, double, double);
    void nearestNeighborLiveUpdateChange(bool value);

    void applyNearestNeighborPressed();
    void hideNearestNeighborPressed();

	// Smoothing
	std::vector<Point3d> smoothPoints(const std::vector<Point3d>& points, Node* rootNode, double radius);
};

#endif
