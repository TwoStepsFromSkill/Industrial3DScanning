#include "RangeQueryWidget.h"

#include <QLabel>
#include <QDoubleSpinBox>
#include <QSpacerItem>
#include <QString>
#include <QPushButton>
#include <QVBoxLayout>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QSlider>
#include <QGroupBox>

#include <iostream>

RangeQueryWidget::RangeQueryWidget(QWidget* parent)
    : QWidget(parent)
{
    m_dummyLayout = new QVBoxLayout(this);
    m_dummyLayout->setContentsMargins(0, 0, 0, 0);
    m_mainWidget = new QGroupBox(this);
    m_mainWidget->setCheckable(true);
    m_dummyLayout->addWidget(m_mainWidget);

    m_mainLayout = new QVBoxLayout(this);
    m_mainWidget->setLayout(m_mainLayout);

    QLabel* positionSectionHeader = new QLabel(QString("Position of query point:"), this);
    positionSectionHeader->setStyleSheet("font-weight: bold");
    m_mainLayout->addWidget(positionSectionHeader);

    QGridLayout* positionLayout = new QGridLayout(this);

    QLabel* xLabel = new QLabel(QString("X"), this);
    xLabel->setStyleSheet("QLabel { color : red; font-weight: bold }");
    QLabel* yLabel = new QLabel(QString("Y"), this);
    yLabel->setStyleSheet("QLabel { color : green; font-weight: bold }");
    QLabel* zLabel = new QLabel(QString("Z"), this);
    zLabel->setStyleSheet("QLabel { color : blue; font-weight: bold }");

    m_xValue = new QSlider(this);
    m_xValue->setMinimum(m_DEFAULT_SLIDER_MIN);
    m_xValue->setMaximum(m_DEFAULT_SLIDER_MAX);
    m_xValue->setTickInterval(1);
    m_xValue->setOrientation(Qt::Horizontal);

    m_yValue = new QSlider(this);
    m_yValue->setMinimum(m_DEFAULT_SLIDER_MIN);
    m_yValue->setMaximum(m_DEFAULT_SLIDER_MAX);
    m_yValue->setTickInterval(1);
    m_yValue->setOrientation(Qt::Horizontal);

    m_zValue = new QSlider(this);
    m_zValue->setMinimum(m_DEFAULT_SLIDER_MIN);
    m_zValue->setMaximum(m_DEFAULT_SLIDER_MAX);
    m_zValue->setTickInterval(1);
    m_zValue->setOrientation(Qt::Horizontal);

    positionLayout->addWidget(xLabel, 0, 0);
    positionLayout->addWidget(m_xValue, 0, 1);

    positionLayout->addWidget(yLabel, 1, 0);
    positionLayout->addWidget(m_yValue, 1, 1);

    positionLayout->addWidget(zLabel, 2, 0);
    positionLayout->addWidget(m_zValue, 2, 1);

    m_positionLabel = new QLabel(QString("x y z"), this);
    positionLayout->addWidget(m_positionLabel, 3, 0, 1, 2);

    m_mainLayout->addLayout(positionLayout);

    QSpacerItem* sectionSpacer = new QSpacerItem(1, 20, QSizePolicy::Fixed, QSizePolicy::Fixed);
    m_mainLayout->addItem(sectionSpacer);

    QLabel* rangeSectionHeader = new QLabel(QString("Range of the query region:"), this);
    rangeSectionHeader->setStyleSheet("font-weight: bold");
    m_mainLayout->addWidget(rangeSectionHeader);

    QGridLayout* rangeLayout = new QGridLayout(this);

    QLabel* dxLabel = new QLabel(QString("X"), this);
    dxLabel->setStyleSheet("QLabel { color : red; font-weight: bold }");
    QLabel* dyLabel = new QLabel(QString("Y"), this);
    dyLabel->setStyleSheet("QLabel { color : green; font-weight: bold }");
    QLabel* dzLabel = new QLabel(QString("Z"), this);
    dzLabel->setStyleSheet("QLabel { color : blue; font-weight: bold }");

    m_dxValue = new QSlider(this);
    m_dxValue->setMinimum(m_DEFAULT_SLIDER_MIN);
    m_dxValue->setMaximum(m_DEFAULT_SLIDER_MAX);
    m_dxValue->setTickInterval(1);
    m_dxValue->setOrientation(Qt::Horizontal);

    m_dyValue = new QSlider(this);
    m_dyValue->setMinimum(m_DEFAULT_SLIDER_MIN);
    m_dyValue->setMaximum(m_DEFAULT_SLIDER_MAX);
    m_dyValue->setTickInterval(1);
    m_dyValue->setOrientation(Qt::Horizontal);

    m_dzValue = new QSlider(this);
    m_dzValue->setMinimum(m_DEFAULT_SLIDER_MIN);
    m_dzValue->setMaximum(m_DEFAULT_SLIDER_MAX);
    m_dzValue->setTickInterval(1);
    m_dzValue->setOrientation(Qt::Horizontal);

    rangeLayout->addWidget(dxLabel, 0, 0);
    rangeLayout->addWidget(m_dxValue, 0, 1);

    rangeLayout->addWidget(dyLabel, 1, 0);
    rangeLayout->addWidget(m_dyValue, 1, 1);

    rangeLayout->addWidget(dzLabel, 2, 0);
    rangeLayout->addWidget(m_dzValue, 2, 1);

    m_rangeLabel = new QLabel(QString("dx dy dz"), this);
    rangeLayout->addWidget(m_rangeLabel, 3, 0, 1, 2);

    m_mainLayout->addLayout(rangeLayout);

    QSpacerItem* buttonSpacer = new QSpacerItem(1, 20, QSizePolicy::Fixed, QSizePolicy::Fixed);
    m_mainLayout->addItem(buttonSpacer);

    QHBoxLayout* buttonLayout = new QHBoxLayout(this);

    m_buttonApply = new QPushButton(QString("Apply"), this);
    m_buttonHide = new QPushButton(QString("Hide"), this);

    buttonLayout->addWidget(m_buttonApply);
    buttonLayout->addWidget(m_buttonHide);

    m_mainLayout->addLayout(buttonLayout);
    m_mainLayout->addStretch(1);

    setLayout(m_dummyLayout);

    connect(m_mainWidget, SIGNAL(clicked(bool)), this, SIGNAL(clickedEnable(bool)));

    connect(m_xValue, SIGNAL(valueChanged(int)), this, SLOT(positionChanged(int)));
    connect(m_yValue, SIGNAL(valueChanged(int)), this, SLOT(positionChanged(int)));
    connect(m_zValue, SIGNAL(valueChanged(int)), this, SLOT(positionChanged(int)));

    connect(m_dxValue, SIGNAL(valueChanged(int)), this, SLOT(rangeChanged(int)));
    connect(m_dyValue, SIGNAL(valueChanged(int)), this, SLOT(rangeChanged(int)));
    connect(m_dzValue, SIGNAL(valueChanged(int)), this, SLOT(rangeChanged(int)));

    connect(m_buttonApply, SIGNAL(pressed()), this, SIGNAL(applyPressed()));
    connect(m_buttonHide, SIGNAL(pressed()), this, SIGNAL(hidePressed()));
}

void RangeQueryWidget::resetValueRange(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax)
{
    m_xValueRange[0] = xMin;
    m_xValueRange[1] = xMax;

    m_yValueRange[0] = yMin;
    m_yValueRange[1] = yMax;

    m_zValueRange[0] = zMin;
    m_zValueRange[1] = zMax;

    m_dxValueRange[0] = 0;
    m_dxValueRange[1] = xMax - xMin;

    m_dyValueRange[0] = 0;
    m_dyValueRange[1] = yMax - yMin;

    m_dzValueRange[0] = 0;
    m_dzValueRange[1] = zMax - zMin;

    updatePositionLabel();
    updateRangeLabel();

    blockSliders(true);

    m_xValue->setValue(0);
    m_yValue->setValue(0);
    m_zValue->setValue(0);

    m_dxValue->setValue(m_INITIAL_RANGE_VALUE);
    m_dyValue->setValue(m_INITIAL_RANGE_VALUE);
    m_dzValue->setValue(m_INITIAL_RANGE_VALUE);

    blockSliders(false);

    emit hidePressed();
    emit centerChanged(getXValue(), getYValue(), getZValue());
    emit extendChanged(getXRange(), getYRange(), getZRange());
}

void RangeQueryWidget::getBox(double* minMax)
{
    minMax[0] = getXValue() - (getXRange() / 2.0);
    minMax[1] = getXValue() + (getXRange() / 2.0);

    minMax[2] = getYValue() - (getYRange() / 2.0);
    minMax[3] = getYValue() + (getYRange() / 2.0);

    minMax[4] = getZValue() - (getZRange() / 2.0);
    minMax[5] = getZValue() + (getZRange() / 2.0);
}

void RangeQueryWidget::positionChanged(int)
{
    updatePositionLabel();
    emit centerChanged(getXValue(), getYValue(), getZValue());
}

void RangeQueryWidget::rangeChanged(int)
{
    updateRangeLabel();
    emit extendChanged(getXRange(), getYRange(), getZRange());
}

double RangeQueryWidget::getXValue() const
{
    return getValue(m_xValue->value(), m_xValueRange);
}

double RangeQueryWidget::getYValue() const
{
    return getValue(m_yValue->value(), m_yValueRange);
}

double RangeQueryWidget::getZValue() const
{
    return getValue(m_zValue->value(), m_zValueRange);
}

double RangeQueryWidget::getXRange() const
{
    return getValue(m_dxValue->value(), m_dxValueRange);
}

double RangeQueryWidget::getYRange() const
{
    return getValue(m_dyValue->value(), m_dyValueRange);
}

double RangeQueryWidget::getZRange() const
{
    return getValue(m_dzValue->value(), m_dzValueRange);
}

double RangeQueryWidget::getValue(int position, const double range[2]) const
{
    double factor = (static_cast<double>(position) - m_DEFAULT_SLIDER_MIN)
                    / (m_DEFAULT_SLIDER_MAX - m_DEFAULT_SLIDER_MIN);
    return (1.0 - factor)*range[0] + factor*range[1];
}

void RangeQueryWidget::blockSliders(bool value)
{
    m_xValue->blockSignals(value);
    m_yValue->blockSignals(value);
    m_zValue->blockSignals(value);

    m_dxValue->blockSignals(value);
    m_dyValue->blockSignals(value);
    m_dzValue->blockSignals(value);
}

void RangeQueryWidget::updatePositionLabel()
{
    m_positionLabel->setText(QString("X: %1 \nY: %2 \nZ: %3").arg(QString::number(getXValue()),
                                                                  QString::number(getYValue()),
                                                                  QString::number(getZValue())));
}

void RangeQueryWidget::updateRangeLabel()
{
    m_rangeLabel->setText(QString("X: %1 \nY: %2 \nZ: %3").arg(QString::number(getXRange()),
                                                               QString::number(getYRange()),
                                                               QString::number(getZRange())));
}