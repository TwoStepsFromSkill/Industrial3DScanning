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

    std::vector<Point3d> m_points;

    Node* m_kdTree;

signals:
    void drawingRangeResultChanged(bool);

private:
    void loadFileXYZ(const char* filename, std::vector<Point3d>& points);
    void updateRangeQueryWidgetData();

private slots:
    void openFile();
    void changeProjection();

    void applyRangePressed();
    void hideRangePressed();
};

#endif
