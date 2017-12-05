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
class BestFitPlaneWidget;

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

    BestFitPlaneWidget* m_bestFitPlaneWidget;

    std::vector<Point3d> m_points;
    Node* m_kdTree;

signals:
    void reloadDrawSettings();

    void drawingMainPointCloudChanged(bool);
    void drawingTemporaryChanged(bool);
    void drawingRangeQueryResultChanged(bool);
    void rangeQueryCenterChange(double, double, double);
    void rangeQueryExtendChange(double, double, double);
    void rangeQueryResultChange(std::vector<Point3d>);

    void drawingNearestNeighborResultChanged(bool);
    void nearestNeighborQueryPointChange(const Point3d&);
    void nearestNeighborResultPointChange(const Point3d&);

    void drawingSmoothedPointsChange(bool);
    void drawingThinnedPointsChange(bool);

    void drawingBestFitPlaneChange(bool);

private:
    void loadFileXYZ(const char* filename, std::vector<Point3d>& points);
    void updateSidebarWidgetData();

    void computeAndVisualizeRangeQuery();
	void computeAndVisualizeNearestNeighbor();

	void computeAndVisualizeSmoothing();
	void computeAndVisualizeThinning();

	std::vector<Point3d> BestFitLine_elke();

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
	void applySmoothing();

    // Thinning
	void applyThinning();

    void smoothTmpPointChanged(int);
    void thinTmpPointChanged(int);

    void applyBestFitPlane();

    std::tuple<Point3d, std::vector<Point3d>, std::vector<Point3d>> bestFitPlane_daniel();
};

#endif
