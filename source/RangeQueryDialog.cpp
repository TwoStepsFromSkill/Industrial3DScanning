#include "RangeQueryDialog.h"

#include <QLabel>
#include <QDoubleSpinBox>
#include <QSpacerItem>
#include <QString>
#include <QPushButton>
#include <QVBoxLayout>
#include <QGridLayout>
#include <QHBoxLayout>

RangeQueryDialog::RangeQueryDialog(QWidget* parent, Qt::WindowFlags f)
    : QDialog(parent, f)
{
    QVBoxLayout* mainLayout = new QVBoxLayout(this);

    QLabel* positionSectionHeader = new QLabel(QString("Position of query point:"), this);
    positionSectionHeader->setStyleSheet("font-weight: bold");
    mainLayout->addWidget(positionSectionHeader);

    QGridLayout* positionLayout = new QGridLayout(this);

    QLabel* xLabel = new QLabel(QString("X Position"), this);
    QLabel* yLabel = new QLabel(QString("Y Position"), this);
    QLabel* zLabel = new QLabel(QString("Z Position"), this);

    m_xValue = new QDoubleSpinBox(this);
    m_xValue->setDecimals(6);
    m_yValue = new QDoubleSpinBox(this);
    m_yValue->setDecimals(6);
    m_zValue = new QDoubleSpinBox(this);
    m_zValue->setDecimals(6);

    m_xValueRangeLabel = new QLabel(QString("(x -- x)"), this);
    m_yValueRangeLabel = new QLabel(QString("(y -- y)"), this);
    m_zValueRangeLabel = new QLabel(QString("(z -- z)"), this);

    positionLayout->addWidget(xLabel, 0, 0);
    positionLayout->addWidget(m_xValue, 0, 1);
    positionLayout->addWidget(m_xValueRangeLabel, 0, 2);

    positionLayout->addWidget(yLabel, 1, 0);
    positionLayout->addWidget(m_yValue, 1, 1);
    positionLayout->addWidget(m_yValueRangeLabel, 1, 2);

    positionLayout->addWidget(zLabel, 2, 0);
    positionLayout->addWidget(m_zValue, 2, 1);
    positionLayout->addWidget(m_zValueRangeLabel, 2, 2);

    mainLayout->addLayout(positionLayout);

    QSpacerItem* spacer = new QSpacerItem(10, 50, QSizePolicy::Minimum, QSizePolicy::Fixed);
    mainLayout->addItem(spacer);

    QLabel* rangeSectionHeader = new QLabel(QString("Range of the query region:"), this);
    rangeSectionHeader->setStyleSheet("font-weight: bold");
    mainLayout->addWidget(rangeSectionHeader);

    QGridLayout* rangeLayout = new QGridLayout(this);

    QLabel* m_dxLabel = new QLabel(QString("X Range"), this);
    QLabel* m_dyLabel = new QLabel(QString("Y Range"), this);
    QLabel* m_dzLabel = new QLabel(QString("Z Range"), this);

    m_dxValue = new QDoubleSpinBox(this);
    m_dxValue->setDecimals(6);
    m_dyValue = new QDoubleSpinBox(this);
    m_dyValue->setDecimals(6);
    m_dzValue = new QDoubleSpinBox(this);
    m_dzValue->setDecimals(6);

    m_dxValueRangeLabel = new QLabel(QString("(x -- x)"), this);
    m_dyValueRangeLabel = new QLabel(QString("(y -- y)"), this);
    m_dzValueRangeLabel = new QLabel(QString("(z -- z)"), this);

    rangeLayout->addWidget(m_dxLabel, 0, 0);
    rangeLayout->addWidget(m_dxValue, 0, 1);
    rangeLayout->addWidget(m_dxValueRangeLabel, 0, 2);

    rangeLayout->addWidget(m_dyLabel, 1, 0);
    rangeLayout->addWidget(m_dyValue, 1, 1);
    rangeLayout->addWidget(m_dyValueRangeLabel, 1, 2);

    rangeLayout->addWidget(m_dzLabel, 2, 0);
    rangeLayout->addWidget(m_dzValue, 2, 1);
    rangeLayout->addWidget(m_dzValueRangeLabel, 2, 2);

    mainLayout->addLayout(rangeLayout);

    QHBoxLayout* buttonLayout = new QHBoxLayout(this);

    m_buttonOK = new QPushButton(QString::fromLatin1("OK"), this);
    m_buttonOK->setDefault(true);
    m_buttonCancel = new QPushButton(QString::fromLatin1("Cancel"), this);

    buttonLayout->addWidget(m_buttonOK);
    buttonLayout->addWidget(m_buttonCancel);

    mainLayout->addLayout(buttonLayout);
    setLayout(mainLayout);

    connect(m_xValue, SIGNAL(valueChanged(double)), this, SLOT(xValueChanged(double)));
    connect(m_yValue, SIGNAL(valueChanged(double)), this, SLOT(yValueChanged(double)));
    connect(m_zValue, SIGNAL(valueChanged(double)), this, SLOT(zValueChanged(double)));

    connect(m_dxValue, SIGNAL(valueChanged(double)), this, SLOT(xRangeChanged(double)));
    connect(m_dyValue, SIGNAL(valueChanged(double)), this, SLOT(yRangeChanged(double)));
    connect(m_dzValue, SIGNAL(valueChanged(double)), this, SLOT(zRangeChanged(double)));

    connect(m_buttonOK, SIGNAL(pressed()), this, SLOT(accept()));
    connect(m_buttonCancel, SIGNAL(pressed()), this, SLOT(reject()));
}

void RangeQueryDialog::setValueRange(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax)
{
    m_xValue->setMinimum(xMin);
    m_xValue->setMaximum(xMax);
    m_xValue->setValue((xMax + xMin) / 2.0);
    m_xValue->setSingleStep((xMax - xMin) / 100);

    m_yValue->setMinimum(yMin);
    m_yValue->setMaximum(yMax);
    m_yValue->setValue((yMax + yMin) / 2.0);
    m_yValue->setSingleStep((xMax - xMin) / 100);

    m_zValue->setMinimum(zMin);
    m_zValue->setMaximum(zMax);
    m_zValue->setValue((zMax + zMin) / 2.0);
    m_zValue->setSingleStep((xMax - xMin) / 100);

    m_dxValue->setMinimum(0);
    m_dxValue->setMaximum(xMax - xMin);
    m_dxValue->setValue((xMax - xMin) / 50.0);
    m_dxValue->setSingleStep((xMax - xMin) / 100);

    m_dyValue->setMinimum(0);
    m_dyValue->setMaximum(yMax - yMin);
    m_dyValue->setValue((yMax - yMin) / 50.0);
    m_dyValue->setSingleStep((xMax - xMin) / 100);

    m_dzValue->setMinimum(0);
    m_dzValue->setMaximum(zMax - zMin);
    m_dzValue->setValue((zMax - zMin) / 50.0);
    m_dzValue->setSingleStep((xMax - xMin) / 100);

    m_xValueRangeLabel->setText(QString("(%1 -- %2)").arg(QString::number(xMin, 'f', 6),
                                                          QString::number(xMax, 'f', 6)));
    m_yValueRangeLabel->setText(QString("(%1 -- %2)").arg(QString::number(yMin, 'f', 6),
                                                          QString::number(yMax, 'f', 6)));
    m_zValueRangeLabel->setText(QString("(%1 -- %2)").arg(QString::number(zMin, 'f', 6),
                                                          QString::number(zMax, 'f', 6)));

    double minRangeX = m_xValue->value() - (m_dxValue->value() / 2.0);
    double maxRangeX = m_xValue->value() - (m_dxValue->value() / 2.0);

    double minRangeY = m_yValue->value() - (m_dyValue->value() / 2.0);
    double maxRangeY = m_yValue->value() - (m_dyValue->value() / 2.0);

    double minRangeZ = m_zValue->value() - (m_dzValue->value() / 2.0);
    double maxRangeZ = m_zValue->value() - (m_dzValue->value() / 2.0);

    m_dxValueRangeLabel->setText(QString("(%1 -- %2)").arg(QString::number(minRangeX, 'f', 6),
                                                           QString::number(maxRangeX, 'f', 6)));
    m_dyValueRangeLabel->setText(QString("(%1 -- %2)").arg(QString::number(minRangeY, 'f', 6),
                                                           QString::number(maxRangeY, 'f', 6)));
    m_dzValueRangeLabel->setText(QString("(%1 -- %2)").arg(QString::number(minRangeZ, 'f', 6),
                                                           QString::number(maxRangeZ, 'f', 6)));
}

void RangeQueryDialog::getPosition(double& x, double& y, double& z)
{
    x = m_xValue->value();
    y = m_yValue->value();
    z = m_zValue->value();
}

void RangeQueryDialog::getRange(double& dx, double& dy, double& dz)
{
    dx = m_dxValue->value();
    dy = m_dyValue->value();
    dz = m_dzValue->value();
}

void RangeQueryDialog::xValueChanged(double newValue)
{
    double newMin = newValue - (m_dxValue->value() / 2.0);
    double newMax = newValue + (m_dxValue->value() / 2.0);

    m_dxValueRangeLabel->setText(QString("(%1 -- %2)").arg(QString::number(newMin, 'f', 6),
                                                           QString::number(newMax, 'f', 6)));
}

void RangeQueryDialog::yValueChanged(double newValue)
{
    double newMin = newValue - (m_dyValue->value() / 2.0);
    double newMax = newValue + (m_dyValue->value() / 2.0);

    m_dyValueRangeLabel->setText(QString("(%1 -- %2)").arg(QString::number(newMin, 'f', 6),
                                                           QString::number(newMax, 'f', 6)));
}

void RangeQueryDialog::zValueChanged(double newValue)
{
    double newMin = newValue - (m_dzValue->value() / 2.0);
    double newMax = newValue + (m_dzValue->value() / 2.0);

    m_dzValueRangeLabel->setText(QString("(%1 -- %2)").arg(QString::number(newMin, 'f', 6),
                                                           QString::number(newMax, 'f', 6)));
}

void RangeQueryDialog::xRangeChanged(double newRange)
{
    double newMin = m_xValue->value() - (newRange / 2.0);
    double newMax = m_xValue->value() - (newRange / 2.0);

    m_dxValueRangeLabel->setText(QString("(%1 -- %2)").arg(QString::number(newMin, 'f', 6),
                                                           QString::number(newMax, 'f', 6)));
}

void RangeQueryDialog::yRangeChanged(double newRange)
{
    double newMin = m_yValue->value() - (newRange / 2.0);
    double newMax = m_yValue->value() - (newRange / 2.0);

    m_dyValueRangeLabel->setText(QString("(%1 -- %2)").arg(QString::number(newMin, 'f', 6),
                                                           QString::number(newMax, 'f', 6)));
}

void RangeQueryDialog::zRangeChanged(double newRange)
{
    double newMin = m_zValue->value() - (newRange / 2.0);
    double newMax = m_zValue->value() - (newRange / 2.0);

    m_dzValueRangeLabel->setText(QString("(%1 -- %2)").arg(QString::number(newMin, 'f', 6),
                                                           QString::number(newMax, 'f', 6)));
}