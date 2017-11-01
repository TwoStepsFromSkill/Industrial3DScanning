#include "MainWindow.h"

#include <QtWidgets/QMenuBar>
#include <QtWidgets/QFileDialog>

#include <QHBoxLayout>
#include <QWidget>
#include <QTabWidget>

#include <fstream>
#include <limits>

#include "RangeQueryWidget.h"

MainWindow::MainWindow()
{
    m_mainLayout = new QHBoxLayout(this);

    m_glWidget = new GLwidget();
    m_glWidget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    m_mainLayout->addWidget(m_glWidget);

    m_tabWidget = new QTabWidget(this);
    m_rangeWidget = new RangeQueryWidget(this);

    m_tabWidget->addTab(m_rangeWidget, QString("Range"));
    m_tabWidget->setTabPosition(QTabWidget::East);
    m_tabWidget->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
    m_tabWidget->setEnabled(false);

    m_mainLayout->addWidget(m_tabWidget);

    QWidget* centralWidget = new QWidget(this);
    centralWidget->setLayout(m_mainLayout);

    setCentralWidget(centralWidget);

    m_fileMenu = menuBar()->addMenu(tr("&File"));
    m_fileMenu->addAction("open",this,SLOT(openFile()));

    m_viewMenu = menuBar()->addMenu(tr("&View"));
    m_viewMenu->addAction("projection", this, SLOT(changeProjection()));

    connect(m_rangeWidget, SIGNAL(clickedEnable(bool)), m_glWidget, SLOT(drawingRangeQueryBoxChange(bool)));

    connect(m_rangeWidget, SIGNAL(centerChanged(double, double, double)),
            m_glWidget, SLOT(rangeQueryCenterChanged(double, double, double)));
    connect(m_rangeWidget, SIGNAL(extendChanged(double, double, double)),
            m_glWidget, SLOT(rangeQueryExtendChanged(double, double, double)));
    connect(m_rangeWidget, SIGNAL(applyPressed()), this, SLOT(applyRangePressed()));
    connect(m_rangeWidget, SIGNAL(hidePressed()), this , SLOT(hideRangePressed()));
    connect(this, SIGNAL(drawingRangeResultChanged(bool)), m_glWidget, SLOT(drawingRangeQueryResultEnabled(bool)));
}

void MainWindow::openFile()
{
    const QString extfilter = ("Pointclouds (*.xyz *.xyzc)");
    QStringList filenames = QFileDialog::getOpenFileNames(this, "Open File", "data", extfilter, 0);

    if (filenames.empty())
        return;

    loadFileXYZ(filenames.front().toLocal8Bit(), m_points);

	Node* tree = KDTree::buildKDTree(m_points.data(), m_points.data() + m_points.size(), 0);
	std::cout << "KDTree was created!" << std::endl;

    m_glWidget->setPoints(m_points);
    m_tabWidget->setEnabled(true);
    updateRangeQueryWidgetData();
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

//Here is the implementation of our file reader
void MainWindow::loadFileXYZ(const char* filename, std::vector<Point3d>& points)
{
    points.clear();
    std::ifstream file(filename);

    if (!file)
    {
        std::cout << "file " << filename << " could not be opened!" << std::endl;
        return; //nothing can be done else -> end function
    }

    std::cout << "reading file: " << filename << std::endl;

    double x = 0, y = 0, z = 0;

    while (file >> x >> y >> z)
    {
        points.emplace_back(x, y, z);
    }

    //dont forget to close to file
    file.close();

    size_t numberOfPoints = points.size();

    std::cout << "reading finished: " << numberOfPoints << " points have be read" << std::endl;
}

void MainWindow::updateRangeQueryWidgetData()
{
    double xMin = std::numeric_limits<double>::max();
    double xMax = std::numeric_limits<double>::lowest();
    double yMin = std::numeric_limits<double>::max();
    double yMax = std::numeric_limits<double>::lowest();
    double zMin = std::numeric_limits<double>::max();
    double zMax = std::numeric_limits<double>::lowest();

    for (const auto& point : m_points)
    {
        xMin = point.x < xMin ? point.x : xMin;
        xMax = point.x > xMax ? point.x : xMax;

        yMin = point.y < yMin ? point.y : yMin;
        yMax = point.y > yMax ? point.y : yMax;

        zMin = point.z < zMin ? point.z : zMin;
        zMax = point.z > zMax ? point.z : zMax;
    }

    m_rangeWidget->resetValueRange(xMin, xMax, yMin, yMax, zMin, zMax);
}

void MainWindow::applyRangePressed()
{
    m_glWidget->setPointsInRange(updateRangeQuery());
    emit drawingRangeResultChanged(true);
}

void MainWindow::hideRangePressed()
{
    emit drawingRangeResultChanged(false);
}

std::vector<Point3d> MainWindow::updateRangeQuery()
{
    std::vector<Point3d> points;

    for (std::size_t i = 0; i < 30; ++i)
        points.push_back(m_points[i]);

    return points;
}