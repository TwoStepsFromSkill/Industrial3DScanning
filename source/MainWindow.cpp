#include "MainWindow.h"

#include <QtWidgets/QMenuBar>
#include <QtWidgets/QFileDialog>

#include <fstream>

MainWindow::MainWindow()
{
    m_glWidget = new GLwidget();
    setCentralWidget(m_glWidget);

    m_fileMenu = menuBar()->addMenu(tr("&File"));
    m_fileMenu->addAction("open",this,SLOT(openFile()));

    m_viewMenu = menuBar()->addMenu(tr("&View"));
    m_viewMenu->addAction("projection", this, SLOT(changeProjection()));
}

void MainWindow::openFile()
{
    const QString extfilter = ("Pointclouds (*.xyz *.xyzc)");
    QStringList filenames = QFileDialog::getOpenFileNames(this, "Open File", "data", extfilter, 0);

    if (filenames.empty()) return;

    std::vector<Point3d> points;
    loadFileXYZ(filenames.front().toLocal8Bit(), points);

	Node* tree = KDTree::buildKDTree(points.data(), points.data() + points.size()-1,0);

	std::cout << tree->leftChild->ptrLastPoint - tree->leftChild->ptrFirstPoint << std::endl;
	std::cout << tree->rightChild->ptrLastPoint - tree->rightChild->ptrFirstPoint << std::endl;

    m_glWidget->setPoints(points);
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
