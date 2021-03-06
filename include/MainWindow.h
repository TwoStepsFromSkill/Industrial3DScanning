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
class BestFitLineWidget;
class BestFitPlaneWidget;
class BestFitSphereWidget;
class NormalWidget;

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

    BestFitLineWidget* m_bestFitLineWidget;
    BestFitPlaneWidget* m_bestFitPlaneWidget;
    BestFitSphereWidget* m_bestFitSphereWidget;

    NormalWidget* m_normalWidget;

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

    void drawingBestFitLineChange(bool);
    void drawingBestFitPlaneChange(bool);
    void drawingBestFitSphereChange(bool);

private:
    void loadFileXYZ(const char* filename, std::vector<Point3d>& points);
    void updateSidebarWidgetData();

    void computeAndVisualizeRangeQuery();
	void computeAndVisualizeNearestNeighbor();

	void computeAndVisualizeSmoothing();
	void computeAndVisualizeThinning();

	std::vector<Point3d> smoothPointsAverage(const std::vector<Point3d>& points, Node* rootNode, double radius);
	/**
	* @brief smoothPointsGaussian is a function that smoothes the points with an gaussian kernel
	in a spherical neighborhood
		* @details  The distance between \f$p` =  \frac{1}{\sum w}\cdot \sum (w_{i}\cdot n_{i})\f$

\f$w_{i} = e^{\frac{-d_{i}}{r}}	\f$

\f$d_{i} = (\sqrt{(p - n_{i})^2})^2\f$.
	* @param points is a vector which contains all given points
	* @param rootNode is a 3dTree which contains all given points in a sorted order
	* @param radius defines the range for the spherical neighboorhood search
	* @return is a vector which contains all smoothed points
	*/
	std::vector<Point3d> smoothPointsGaussian(const std::vector<Point3d>& points, Node* rootNode, double radius);

    double stdDeviation(const std::vector<double>& values) const;

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

    void applyBestFitLine();
    void applyBestFitPlane();
    void applyBestFitSphere();

    void applyVertexNormal();

    std::tuple<Point3d, std::vector<Point3d>, std::vector<Point3d>,
               std::vector<Point3d>> bestFitPlane_daniel();

	std::tuple<Point3d, std::vector<Point3d>, std::vector<Point3d>,
		std::vector<Point3d>> MainWindow::bestFitPlaneForSomePoints(std::vector<Point3d> points);
	/**
	* @brief bestFitSphere_elke is a function that computes the best fit sphere for a given point cloud m_points
	* @details  The parameter for the function of the sphere were iteratively computed by a singular value decomposition
	and a linear equation system
	* @return is a vector which contains all parameters of the function of the sphere
	*/
	std::vector<double> bestFitSphere_elke();

    std::vector<Point3d> computeVisualSphere(const Point3d& center, double r);
	
	std::vector<Point3d> calculateAllNormals(const double radius);
	
	Point3d calculateNormalVector(std::vector<Point3d> plane);
};

#endif
