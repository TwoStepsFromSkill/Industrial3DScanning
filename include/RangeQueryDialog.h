#ifndef RANGE_QUERY_DIALOG_H
#define RANGE_QUERY_DIALOG_H

#include <QDialog>

class QLabel;
class QDoubleSpinBox;
class QPushButton;

class RangeQueryDialog : public QDialog
{
    Q_OBJECT
public:
    RangeQueryDialog(QWidget* parent = nullptr, Qt::WindowFlags f = Qt::WindowFlags());

    void setValueRange(double xMin, double xMax,
                       double yMin, double yMax,
                       double zMin, double zMax);

    void getPosition(double& x, double& y, double& z);
    void getRange(double& dx, double& dy, double& dz);

private:
    QDoubleSpinBox* m_xValue;
    QDoubleSpinBox* m_yValue;
    QDoubleSpinBox* m_zValue;

    QLabel* m_xValueRangeLabel;
    QLabel* m_yValueRangeLabel;
    QLabel* m_zValueRangeLabel;

    QDoubleSpinBox* m_dxValue;
    QDoubleSpinBox* m_dyValue;
    QDoubleSpinBox* m_dzValue;

    QLabel* m_dxValueRangeLabel;
    QLabel* m_dyValueRangeLabel;
    QLabel* m_dzValueRangeLabel;

    QPushButton* m_buttonOK;
    QPushButton* m_buttonCancel;

private slots:
    void xValueChanged(double newValue);
    void yValueChanged(double newValue);
    void zValueChanged(double newValue);

    void xRangeChanged(double newRange);
    void yRangeChanged(double newRange);
    void zRangeChanged(double newRange);
};

#endif // RANGE_QUERY_DIALOG_H
