#include "MainWindow.h"

#include <QtWidgets/QMenuBar>
#include <QtWidgets/QFileDialog>

#include <fstream>
#include <limits>

MainWindow::MainWindow()
{
    m_glWidget = new GLwidget();
    setCentralWidget(m_glWidget);

    m_fileMenu = menuBar()->addMenu(tr("&File"));
    m_fileMenu->addAction("open",this,SLOT(openFile()));

    m_viewMenu = menuBar()->addMenu(tr("&View"));
    m_viewMenu->addAction("projection", this, SLOT(changeProjection()));

    m_rangeMenu = menuBar()->addMenu(tr("&Query"));
    m_rangeMenu->addAction("range query", this, SLOT(setRangeQuery()));

    m_rangeDialog = new RangeQueryDialog(this);
}

void MainWindow::openFile()
{
    const QString extfilter = ("Pointclouds (*.xyz *.xyzc)");
    QStringList filenames = QFileDialog::getOpenFileNames(this, "Open File", "data", extfilter, 0);

    if (filenames.empty()) return;

    loadFileXYZ(filenames.front().toLocal8Bit(), m_points);

	Node* tree = KDTree::buildKDTree(m_points.data(), m_points.data() + m_points.size(),0);

	std::cout << "KDTree was created!" << std::endl;

    m_glWidget->setPoints(m_points);
}

void MainWindow::changeProjection()
{
    printf("change projection\n");

    m_glWidget->makeCurrent();

    if (m_glWidget->camera().usesPerspectiveProjection())
        m_glWidget->camera().usePerspectiveProjection(false);
    else
        m_glWidget->camera().usePerspectiveProjection(true);

    m_glWidget->update();
}

void MainWindow::setRangeQuery()
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

    m_rangeDialog->setValueRange(xMin, xMax, yMin, yMax, zMin, zMax);

    if (m_rangeDialog->exec())
    {
        double xPos, yPos, zPos;
        double dx, dy, dz;

        m_rangeDialog->getPosition(xPos, yPos, zPos);
        m_rangeDialog->getRange(dx, dy, dz);

        std::cout << "Position: (" << xPos << ", " << yPos << ", " << zPos << ")\n";
        std::cout << "Range: " << dx << ", " << dy << ", " << dz << "\n\n";
    }
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

    double x = 0;
    double y = 0;
    double z = 0;

    while (file >> x >> y >> z)
    {
        points.emplace_back(x, y, z);
    }

    //dont forget to close to file
    file.close();

    size_t numberOfPoints = points.size();

    std::cout << "reading finished: " << numberOfPoints << " points have be read" << std::endl;
}
