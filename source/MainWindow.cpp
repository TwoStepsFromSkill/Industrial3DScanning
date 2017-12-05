#include "MainWindow.h"

#include <QtWidgets/QMenuBar>
#include <QtWidgets/QFileDialog>

#include <QHBoxLayout>
#include <QWidget>
#include <QTabWidget>

#include <fstream>
#include <limits>
#include <chrono>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <array>
#include <utility>

#include "BaseTabWidget.h"
#include "RangeQueryWidget.h"
#include "NearestNeighborWidget.h"
#include "SmoothingWidget.h"
#include "ThinningWidget.h"
#include "KDTree.h"
#include "BestFitPlaneWidget.h"
#include "Matrix.h"
#include "SVD.h"

using duration_micro = std::chrono::duration<double, std::micro>;
using duration_milli = std::chrono::duration<double, std::milli>;

MainWindow::MainWindow()
    : m_liveUpdateRangeQuery(false)
    , m_liveUpdateNearestNeighbor(false)
    , m_points()
    , m_kdTree(nullptr)
{
    m_mainLayout = new QHBoxLayout();

    m_glWidget = new GLwidget();
    m_glWidget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    this->installEventFilter(m_glWidget);

    m_mainLayout->addWidget(m_glWidget);

    m_tabWidget = new QTabWidget(this);
    m_tabWidget->setTabPosition(QTabWidget::East);
    m_tabWidget->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
    m_tabWidget->setEnabled(false);

    m_rangeWidget = new RangeQueryWidget(this);
    m_tabWidget->addTab(m_rangeWidget, QString("Range Query"));

    m_nearestWidget = new NearestNeighborWidget(this);
    m_tabWidget->addTab(m_nearestWidget, QString("Nearest Neighboor"));

    m_smoothingWidget = new SmoothingWidget(this);
    m_tabWidget->addTab(m_smoothingWidget, QString("Smoothing"));

    m_thinningWidget = new ThinningWidget(this);
    m_tabWidget->addTab(m_thinningWidget, QString("Thinning"));

    m_bestFitPlaneWidget = new BestFitPlaneWidget(this);
    m_tabWidget->addTab(m_bestFitPlaneWidget, QString("Best Fit Plane"));

    m_mainLayout->addWidget(m_tabWidget);

    QWidget* centralWidget = new QWidget(this);
    centralWidget->setLayout(m_mainLayout);

    setCentralWidget(centralWidget);

    m_fileMenu = menuBar()->addMenu(tr("&File"));
    m_fileMenu->addAction("Open XYZ-File", this, SLOT(openFile()));

    m_viewMenu = menuBar()->addMenu(tr("&View"));
    m_viewMenu->addAction("Toggle Projection", this, SLOT(changeProjection()));
    m_viewMenu->addAction("Reload draw settings", this, SIGNAL(reloadDrawSettings()));

    connect(this, SIGNAL(reloadDrawSettings()), m_glWidget, SLOT(reloadDrawSettings()));

    // TabWidget section
    connect(m_tabWidget, SIGNAL(currentChanged(int)), this, SLOT(tabSwitched(int)));

    // RangeQueryWidget section
    connect(m_rangeWidget, SIGNAL(widgetEnabled(bool)),
            m_glWidget, SLOT(drawingRangeQueryBoxChange(bool)));

    connect(m_rangeWidget, SIGNAL(centerChanged(double, double, double)),
            this, SLOT(rangeQueryCenterChanged(double, double, double)));
    connect(m_rangeWidget, SIGNAL(extendChanged(double, double, double)),
            this, SLOT(rangeQueryExtendChanged(double, double, double)));

    connect(this, SIGNAL(rangeQueryCenterChange(double, double, double)),
            m_glWidget, SLOT(rangeQueryCenterChanged(double, double, double)));
    connect(this, SIGNAL(rangeQueryExtendChange(double, double, double)),
            m_glWidget, SLOT(rangeQueryExtendChanged(double, double, double)));

    connect(m_rangeWidget, SIGNAL(liveUpdateChanged(bool)),
            this, SLOT(rangeQueryLiveUpdateChange(bool)));

    connect(m_rangeWidget, SIGNAL(applyPressed()), this, SLOT(applyRangeQueryPressed()));
    connect(m_rangeWidget, SIGNAL(hidePressed()), this , SLOT(hideRangeQueryPressed()));

    connect(this, SIGNAL(drawingRangeQueryResultChanged(bool)),
            m_glWidget, SLOT(drawingRangeQueryResultEnabled(bool)));
    connect(this, SIGNAL(rangeQueryResultChange(std::vector<Point3d>)),
            m_glWidget, SLOT(rangeQueryResultChanged(std::vector<Point3d>)));

    // NearestNeighborWidget section
    connect(m_nearestWidget, SIGNAL(widgetEnabled(bool)),
            m_glWidget, SLOT(drawingNearestNeighborQueryPointChanged(bool)));

    connect(m_nearestWidget, SIGNAL(positionChange(double, double, double)),
            this, SLOT(nearestNeighborQueryPointChanged(double, double, double)));
    connect(this, SIGNAL(nearestNeighborQueryPointChange(const Point3d&)),
            m_glWidget, SLOT(nearestNeighborQueryPointChanged(const Point3d&)));

    connect(m_nearestWidget, SIGNAL(liveUpdateChanged(bool)),
            this, SLOT(nearestNeighborLiveUpdateChange(bool)));

    connect(m_nearestWidget, SIGNAL(applyPressed()), this, SLOT(applyNearestNeighborPressed()));
    connect(m_nearestWidget, SIGNAL(hidePressed()), this , SLOT(hideNearestNeighborPressed()));

    connect(this, SIGNAL(drawingNearestNeighborResultChanged(bool)),
            m_glWidget, SLOT(drawingNearestNeighborResultPointChanged(bool)));
    connect(this, SIGNAL(nearestNeighborResultPointChange(const Point3d&)),
            m_glWidget, SLOT(nearestNeighborResultPointChanged(const Point3d&)));

    m_rangeWidget->activate();

	// Smoothing section
	connect(m_smoothingWidget, SIGNAL(applyPressed()), this, SLOT(applySmoothing()));
    connect(this, SIGNAL(drawingSmoothedPointsChange(bool)), m_glWidget, SLOT(drawingSmoothedPointsChanged(bool)));
    connect(m_smoothingWidget, SIGNAL(tempPointChanged(int)), this, SLOT(smoothTmpPointChanged(int)));

    // Thinning section
    connect(m_thinningWidget, SIGNAL(applyPressed()), this, SLOT(applyThinning()));
    connect(this, SIGNAL(drawingThinnedPointsChange(bool)), m_glWidget, SLOT(drawingThinnedPointsChanged(bool)));
    connect(m_thinningWidget, SIGNAL(tempPointChanged(int)), this, SLOT(thinTmpPointChanged(int)));

    connect(this, SIGNAL(drawingMainPointCloudChanged(bool)), m_glWidget, SLOT(drawingMainPointCloudChange(bool)));
    connect(this, SIGNAL(drawingTemporaryChanged(bool)), m_glWidget, SLOT(drawingTemporaryChange(bool)));

    // Best Fit Plane
    connect(m_bestFitPlaneWidget, SIGNAL(applyPressed()), this, SLOT(applyBestFitPlane()));
    connect(this, SIGNAL(drawingBestFitPlaneChange(bool)), m_glWidget, SLOT(drawingBestFitPlaneChanged(bool)));
}

void MainWindow::openFile()
{
    const QString extfilter = ("Pointclouds (*.xyz *.xyzc)");
    QStringList filenames = QFileDialog::getOpenFileNames(this, "Open File", "data", extfilter, 0);

    if (filenames.empty())
        return;

    auto startTime = std::chrono::system_clock::now();
        loadFileXYZ(filenames.front().toLocal8Bit(), m_points);
    duration_milli elapsed = std::chrono::system_clock::now() - startTime;
    std::cerr << "Loaded " << m_points.size() << " poinst! Took [" << elapsed.count() << "ms]\n";

    if (m_kdTree)
    {
        delete m_kdTree;
    }

    startTime = std::chrono::system_clock::now();
        m_kdTree = KDTree::buildKDTree(m_points.data(), m_points.data() + m_points.size(), 0);
    elapsed = std::chrono::system_clock::now() - startTime;
    std::cerr << "KDTree created! Took [" << elapsed.count() << "ms]\n";

    m_glWidget->setPoints(m_points);
    m_tabWidget->setEnabled(true);

    updateSidebarWidgetData();
    m_glWidget->drawingRangeQueryBoxChange(true);

}

void MainWindow::changeProjection()
{
    m_glWidget->makeCurrent();

    if (m_glWidget->camera().usesPerspectiveProjection())
        m_glWidget->camera().usePerspectiveProjection(false);
    else
        m_glWidget->camera().usePerspectiveProjection(true);

    m_glWidget->update();
}

void MainWindow::tabSwitched(int index)
{
    for (int i = 0; i < m_tabWidget->count(); ++i)
    {
        if (i != index)
           static_cast<BaseTabWidget*>(m_tabWidget->widget(i))->deactivate();
        else
            static_cast<BaseTabWidget*>(m_tabWidget->widget(i))->activate();
    }

    emit drawingTemporaryChanged(false);
    emit drawingSmoothedPointsChange(false);
    emit drawingThinnedPointsChange(false);
}

void MainWindow::loadFileXYZ(const char* filename, std::vector<Point3d>& points)
{
    points.clear();
    std::ifstream file(filename);

    if (!file)
    {
        std::cerr << "file " << filename << " could not be opened!" << std::endl;
        return;
    }

    double x = 0, y = 0, z = 0;

    while (file >> x >> y >> z)
    {
        points.emplace_back(x, y, 0.0);
    }

    file.close();
}

void MainWindow::updateSidebarWidgetData()
{
    double xMin = std::numeric_limits<double>::max();
    double xMax = std::numeric_limits<double>::lowest();
    double yMin = std::numeric_limits<double>::max();
    double yMax = std::numeric_limits<double>::lowest();
    double zMin = std::numeric_limits<double>::max();
    double zMax = std::numeric_limits<double>::lowest();

    for (const auto& point : m_points)
    {
        xMin = point[0] < xMin ? point[0] : xMin;
        xMax = point[0] > xMax ? point[0] : xMax;

        yMin = point[1] < yMin ? point[1] : yMin;
        yMax = point[1] > yMax ? point[1] : yMax;

        zMin = point[2] < zMin ? point[2] : zMin;
        zMax = point[2] > zMax ? point[2] : zMax;
    }

    m_rangeWidget->resetValueRange(xMin, xMax, yMin, yMax, zMin, zMax);
    m_nearestWidget->resetValueRange(xMin, xMax, yMin, yMax, zMin, zMax);
    m_smoothingWidget->resetValueRange(xMin, xMax, yMin, yMax, zMin, zMax);
    m_smoothingWidget->setNumberOfPoints(m_points.size());
    m_thinningWidget->resetValueRange(xMin, xMax, yMin, yMax, zMin, zMax);
    m_thinningWidget->setNumberOfPoints(m_points.size());
}

void MainWindow::rangeQueryCenterChanged(double x, double y, double z)
{
    emit rangeQueryCenterChange(x, y, z);

    if (m_liveUpdateRangeQuery)
    {
        computeAndVisualizeRangeQuery();
    }
}

void MainWindow::rangeQueryExtendChanged(double dx, double dy, double dz)
{
    emit rangeQueryExtendChange(dx, dy, dz);

    if (m_liveUpdateRangeQuery)
    {
        computeAndVisualizeRangeQuery();
    }
}

void MainWindow::rangeQueryLiveUpdateChange(bool value)
{
    m_liveUpdateRangeQuery = value;
}

void MainWindow::applyRangeQueryPressed()
{
    if (m_kdTree)
    {
        computeAndVisualizeRangeQuery();
    }
}

void MainWindow::hideRangeQueryPressed()
{
    emit drawingRangeQueryResultChanged(false);
}

void MainWindow::computeAndVisualizeRangeQuery()
{
    double range[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    m_rangeWidget->getBox(&range[0]);

    auto startTime = std::chrono::system_clock::now();
    auto points = queryRange(m_kdTree, range);
    duration_micro elapsed = std::chrono::system_clock::now() - startTime;
    std::cerr << "Queried KDTree! Found " << points.size() << " poinst! Took [" << elapsed.count() << "µs]\n";

    emit rangeQueryResultChange(points);
    emit drawingRangeQueryResultChanged(true);
}

void MainWindow::nearestNeighborQueryPointChanged(double x, double y, double z)
{
    emit nearestNeighborQueryPointChange(Point3d(x, y, z));

    if (m_liveUpdateNearestNeighbor)
    {
        computeAndVisualizeNearestNeighbor();
    }
}

void MainWindow::nearestNeighborLiveUpdateChange(bool value)
{
    m_liveUpdateNearestNeighbor = value;
}

void MainWindow::applyNearestNeighborPressed()
{
    if (m_kdTree)
    {
        computeAndVisualizeNearestNeighbor();
    }
}

void MainWindow::hideNearestNeighborPressed()
{
    emit drawingNearestNeighborResultChanged(false);
}

void MainWindow::computeAndVisualizeNearestNeighbor()
{
    double xyz[3];
    m_nearestWidget->getQueryPoint(xyz);

    auto startTime = std::chrono::system_clock::now();
        Point3d res = nearestNeighbor_daniel(m_kdTree, Point3d(xyz[0], xyz[1], xyz[2]));
    duration_micro elapsed = std::chrono::system_clock::now() - startTime;
    std::cerr << "Found NearestNeighbor! Took [" << elapsed.count() << "µs]\n";

    emit nearestNeighborResultPointChange(res);
    emit drawingNearestNeighborResultChanged(true);
}

void MainWindow::applySmoothing()
{
	if (m_kdTree)
	{
		computeAndVisualizeSmoothing();
	}
}

void MainWindow::applyThinning()
{
	if (m_kdTree)
	{
		computeAndVisualizeThinning();
	}
}

void MainWindow::smoothTmpPointChanged(int i)
{
    const Point3d& pt = m_points[i];
    double radius = 0.0;
    m_smoothingWidget->getRadius(&radius);

    std::vector<Point3d> neighborPoints = queryRadius(m_kdTree, radius, pt);

    m_glWidget->setTempPoint(pt);
    m_glWidget->setTempRadiusPoints(neighborPoints);
    emit drawingTemporaryChanged(true);
}

void MainWindow::thinTmpPointChanged(int i)
{
    const Point3d& pt = m_points[i];
    double radius = 0.0;
    m_thinningWidget->getRadius(&radius);

    std::vector<Point3d> neighborPoints = queryRadius(m_kdTree, radius, pt);

    m_glWidget->setTempPoint(pt);
    m_glWidget->setTempRadiusPoints(neighborPoints);
    emit drawingTemporaryChanged(true);
}

void MainWindow::computeAndVisualizeSmoothing()
{
	double radius;
	m_smoothingWidget->getRadius(&radius);

    std::vector<Point3d> smoothedPoints;

    auto startTime = std::chrono::system_clock::now();
    if (m_smoothingWidget->useGaussSmoothing())
    {
        smoothedPoints = smoothPointsGaussian(m_points, m_kdTree, radius);
    }
    else
    {
        smoothedPoints = smoothPointsAverage(m_points, m_kdTree, radius);
    }
    duration_milli elapsed = std::chrono::system_clock::now() - startTime;
    std::cerr << "Finished Smoothing! Took [" << elapsed.count() / 1000.0 << "s]\n";

    std::vector<double> distances(m_points.size());

#pragma omp parallel for
    for (int i = 0; i < m_points.size(); ++i)
    {
        distances[i] = sqDistance3d(m_points[i], smoothedPoints[i]);
    }

    m_glWidget->setSmoothedPoints(smoothedPoints);
    m_glWidget->setPointDistances(distances);

    emit drawingTemporaryChanged(false);
    emit drawingMainPointCloudChanged(false);
    emit drawingSmoothedPointsChange(true);
}

void MainWindow::computeAndVisualizeThinning()
{
    double radius;
	m_thinningWidget->getRadius(&radius);

#pragma omp parallel for
    for (int i = 0; i < m_points.size(); ++i)
    {
        m_points[i].flag_ignore = false;
    }

    std::vector<Point3d> thinnedPoints;

    std::cerr << "Started thinning...\n";
    auto startTime = std::chrono::system_clock::now();
        homogeneousThinning(m_kdTree, m_kdTree, radius, thinnedPoints);
    duration_milli elapsed = std::chrono::system_clock::now() - startTime;
    std::cerr << "Finished Thinning! Took [" << elapsed.count() / 1000.0 << "s]\n";

    m_glWidget->setThinnedPoints(thinnedPoints);
    emit drawingTemporaryChanged(false);
    emit drawingMainPointCloudChanged(false);
    emit drawingThinnedPointsChange(true);
}


void MainWindow::applyBestFitPlane()
{
    auto startTime = std::chrono::system_clock::now();
        auto planeParts = bestFitPlane_daniel();
    
	duration_milli elapsed = std::chrono::system_clock::now() - startTime;
    std::cerr << "Computed BFP! Took [" << elapsed.count() << "ms]\n";
	
	auto startTimeBFL = std::chrono::system_clock::now();
	auto linePoints = BestFitLine_elke();
	elapsed = std::chrono::system_clock::now() - startTimeBFL;
	std::cerr << "Computed BFL! Took [" << elapsed.count() << "ms]\n";

    const Point3d& C = std::get<0>(planeParts);
    const Point3d& EV = std::get<2>(planeParts)[2];
    std::vector<double> distances(m_points.size());

#pragma omp parallel for
    for (int i = 0; i < m_points.size(); ++i)
    {
        distances[i] = dotProduct(EV, m_points[i] - C);
    }

    m_glWidget->setBFPCorners(std::get<1>(planeParts));
    m_glWidget->setPointDistances(distances);
    m_glWidget->drawingMainPointWithColorArray(true);

	m_glWidget->setBFLPoints(linePoints);

    emit drawingBestFitPlaneChange(true);
}

std::tuple<Point3d, std::vector<Point3d>, std::vector<Point3d>> MainWindow::bestFitPlane_daniel()
{
    // Computer center (mean)
    double centerX = 0;
    double centerY = 0;
    double centerZ = 0;

#pragma omp parallel for reduction(+:centerX,centerY,centerZ)
    for (int i = 0; i < m_points.size(); ++i)
    {
        centerX += m_points[i][0];
        centerY += m_points[i][1];
        centerZ += m_points[i][2];
    }

    centerX /= m_points.size();
    centerY /= m_points.size();
    centerZ /= m_points.size();

    // Compute covariance matrix
    double Cxx = 0; double Cxy = 0; double Cxz = 0;
    double Cyy = 0; double Cyz = 0;
    double Czz = 0;

    std::size_t n = m_points.size() - 1;

#pragma omp parallel for reduction(+:centerX,centerY,centerZ)
    for (int i = 0; i < m_points.size(); ++i)
    {
        Cxx += (m_points[i][0] - centerX)*(m_points[i][0] - centerX);
        Cxy += (m_points[i][0] - centerX)*(m_points[i][1] - centerY);
        Cxz += (m_points[i][0] - centerX)*(m_points[i][2] - centerZ);

        Cyy += (m_points[i][1] - centerY)*(m_points[i][1] - centerY);
        Cyz += (m_points[i][1] - centerY)*(m_points[i][2] - centerZ);

        Czz += (m_points[i][2] - centerZ)*(m_points[i][2] - centerZ);
    }

    Matrix cov(3,3);

    cov(0,0) = Cxx / n;
    cov(0,1) = Cxy / n;
    cov(0,2) = Cxz / n;

    cov(1,0) = cov(0,1);
    cov(1,1) = Cyy / n;
    cov(1,2) = Cyz / n;

    cov(2,0) = cov(0,2);
    cov(2,1) = cov(1,2);
    cov(2,2) = Czz / n;

    SVD::computeSymmetricEigenvectors(cov);

    Point3d EV0(cov(0,0), cov(1,0), cov(2,0));
    Point3d EV1(cov(0,1), cov(1,1), cov(2,1));
    Point3d EV2(cov(0,2), cov(1,2), cov(2,2));

    normalizeVector(EV0);
    normalizeVector(EV1);
    normalizeVector(EV2);

    std::vector<Point3d> corners;
    Point3d C(centerX, centerY, centerZ);

    double maxDistEV0 = std::numeric_limits<double>::lowest();
    double minDistEV0 = std::numeric_limits<double>::max();

    double maxDistEV1 = std::numeric_limits<double>::lowest();
    double minDistEV1 = std::numeric_limits<double>::max();

    for (std::size_t i = 0; i < m_points.size(); ++i)
    {
        double dist = dotProduct(EV0, m_points[i] - C);
        maxDistEV0 = dist > maxDistEV0 ? dist :  maxDistEV0;
        minDistEV0 = dist < minDistEV0 ? dist :  minDistEV0;

        dist = dotProduct(EV1, m_points[i] - C);
        maxDistEV1 = dist > maxDistEV1 ? dist :  maxDistEV1;
        minDistEV1 = dist < minDistEV1 ? dist :  minDistEV1;
    }

    corners.push_back(C + EV0*maxDistEV0 + EV1*maxDistEV1);
    corners.push_back(C + EV0*maxDistEV0 + EV1*minDistEV1);
    corners.push_back(C + EV0*minDistEV0 + EV1*minDistEV1);
    corners.push_back(C + EV0*minDistEV0 + EV1*maxDistEV1);

    std::vector<Point3d> evs;
    evs.push_back(EV0);
    evs.push_back(EV1);
    evs.push_back(EV2);

    return std::make_tuple(C, corners, evs);
}

std::vector<Point3d> MainWindow::BestFitLine_elke()
{
	// Computer center (mean)
	double centerX = 0;
	double centerY = 0;
	double centerZ = 0;

#pragma omp parallel for reduction(+:centerX,centerY,centerZ)
	for (int i = 0; i < m_points.size(); ++i)
	{
		centerX += m_points[i][0];
		centerY += m_points[i][1];
		centerZ += m_points[i][2];
	}

	centerX /= m_points.size();
	centerY /= m_points.size();
	centerZ /= m_points.size();

	// Compute covariance matrix
	double Cxx = 0; double Cxy = 0; double Cxz = 0;
	double Cyy = 0; double Cyz = 0;
	double Czz = 0;

	std::size_t n = m_points.size() - 1;

#pragma omp parallel for reduction(+:centerX,centerY,centerZ)
	for (int i = 0; i < m_points.size(); ++i)
	{
		Cxx += (m_points[i][0] - centerX)*(m_points[i][0] - centerX);
		Cxy += (m_points[i][0] - centerX)*(m_points[i][1] - centerY);
		Cxz += (m_points[i][0] - centerX)*(m_points[i][2] - centerZ);

		Cyy += (m_points[i][1] - centerY)*(m_points[i][1] - centerY);
		Cyz += (m_points[i][1] - centerY)*(m_points[i][2] - centerZ);

		Czz += (m_points[i][2] - centerZ)*(m_points[i][2] - centerZ);
	}

	Matrix cov(3, 3);

	cov(0, 0) = Cxx / n;
	cov(0, 1) = Cxy / n;
	cov(0, 2) = Cxz / n;

	cov(1, 0) = cov(0, 1);
	cov(1, 1) = Cyy / n;
	cov(1, 2) = Cyz / n;

	cov(2, 0) = cov(0, 2);
	cov(2, 1) = cov(1, 2);
	cov(2, 2) = Czz / n;

	SVD::computeSymmetricEigenvectors(cov);

	Point3d EV0(cov(0, 0), cov(1, 0), cov(2, 0));

	normalizeVector(EV0);

	Point3d C(centerX, centerY, centerZ);

	double maxDistEV0 = std::numeric_limits<double>::lowest();
	double minDistEV0 = std::numeric_limits<double>::max();

	for (std::size_t i = 0; i < m_points.size(); ++i)
	{
		double dist = dotProduct(EV0, m_points[i] - C);
		maxDistEV0 = dist > maxDistEV0 ? dist : maxDistEV0;
		minDistEV0 = dist < minDistEV0 ? dist : minDistEV0;
	}
	std::vector<Point3d> line;
	line.push_back(C + EV0*minDistEV0);
	line.push_back(C + EV0*maxDistEV0);

	return line;
}