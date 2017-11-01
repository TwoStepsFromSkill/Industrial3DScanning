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

	//----TEST---AREA----
	//get points from file
	openFile();

	//create K-d-Tree from points
	Node* testTree = KDTree::buildKDTree(m_points.data(), m_points.data() + m_points.size(), 0);
	
	/*std::cout << "testTree.ptrFirstPoint.x = " + std::to_string(testTree->ptrFirstPoint->x) + "\n";
	std::cout << "testTree.ptrFirstPoint.y = " + std::to_string(testTree->ptrFirstPoint->y) + "\n";
	std::cout << "testTree.ptrFirstPoint.z = " + std::to_string(testTree->ptrFirstPoint->z) + "\n\n";

	std::cout << "testTree.ptrLastPoint.x = " + std::to_string(testTree->ptrLastPoint->x) + "\n";
	std::cout << "testTree.ptrLastPoint.y = " + std::to_string(testTree->ptrLastPoint->y) + "\n";
	std::cout << "testTree.ptrLastPoint.z = " + std::to_string(testTree->ptrLastPoint->z) + "\n\n";*/

	//get points in specified range
	Point3d rangeMin = Point3d(0.0, 0.0, 0.0);
	Point3d rangeMax = Point3d(5.0, 5.0, 5.0);
	std::vector<Point3d> testQuery = getPointsInRange(testTree, &rangeMin, &rangeMax);
}

std::vector<Point3d> MainWindow::getPointsInRange(Node* KDTree, Point3d* rMin, Point3d* rMax)
{
	//create vector to save query results
	std::vector<Point3d> resVector;

	//search KDTree for points in range and save them in resVector
	searchKDTree(KDTree, 0, rMin, rMax, &resVector);

	//return vector with query results
	return resVector;
}

void MainWindow::searchKDTree(Node* curNode, unsigned int depth, Point3d* rMin, Point3d* rMax, std::vector<Point3d>* queryVector)
{
	if (curNode) std::cout << "\ncurNode = true\n";
	if (curNode->median) std::cout << "curNode.median = true, Value = " << curNode->median << "\n";
	if (curNode->leftChild) std::cout << "curNode.leftChild = true\n";
	if (curNode->rightChild) std::cout << "curNode.rightChild = true\n";
	if (curNode->ptrFirstPoint) std::cout << "curNode.ptrFirstPoint = true, " << Point3dToString(curNode->ptrFirstPoint) << "\n";
	if (curNode->ptrLastPoint) std::cout << "curNode.ptrLastPoint = true, " << Point3dToString(curNode->ptrLastPoint) << "\n";

	//calculate current dimension via Modulo-Operation
	unsigned int currentDimension = depth % 3;

	//compare median value depending on the current tree depth
	if (currentDimension == 0)
	{
		//std::cout << "Depth: " << depth << ", X-Median = " << curNode->median << "\n";
		if (rMin->x < curNode->median)
		{
			std::cout << "-> going left" << "\n";
			searchKDTree(curNode->leftChild, depth + 1, rMin, rMax, queryVector);
		}
		else if (rMax->x >= curNode->median)
		{
			std::cout << "-> going right" << "\n";
			searchKDTree(curNode->rightChild, depth + 1, rMin, rMax, queryVector);
		}
	}
	else if (currentDimension == 1)
	{
		//std::cout << "Depth: " << depth << ", Y-Median = " << curNode->median << "\n";
		if (rMin->y < curNode->median)
		{
			std::cout << "-> going left" << "\n";
			searchKDTree(curNode->leftChild, depth + 1, rMin, rMax, queryVector);
		}
		else if (rMax->y >= curNode->median)
		{
			std::cout << "-> going right" << "\n";
			searchKDTree(curNode->rightChild, depth + 1, rMin, rMax, queryVector);
		}
	}
	else
	{
		//std::cout << "Depth: " << depth << ", Z-Median = " << curNode->median << "\n";
		if (rMin->z < curNode->median)
		{
			std::cout << "-> going left" << "\n";
			searchKDTree(curNode->leftChild, depth + 1, rMin, rMax, queryVector);
		}
		else if (rMax->z >= curNode->median)
		{
			std::cout << "-> going right" << "\n";
			searchKDTree(curNode->rightChild, depth + 1, rMin, rMax, queryVector);
		}
	}

	//add point to resVector
	std::cout << "Query: " << Point3dToString(curNode->ptrFirstPoint) << "\n";
	queryVector->push_back(*(curNode->ptrFirstPoint));
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
