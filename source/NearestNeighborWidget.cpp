#include "NearestNeighborWidget.h"

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
#include <QCheckBox>

#include <iostream>

NearestNeighborWidget::NearestNeighborWidget(QWidget* parent)
    : BaseTabWidget(parent)
{
    m_dummyLayout = new QVBoxLayout();
    m_dummyLayout->setContentsMargins(0, 0, 0, 0);
    m_mainWidget = new QGroupBox(this);
    m_mainWidget->setCheckable(true);
    m_mainWidget->setTitle(QString("Nearest Neighboor"));
    m_dummyLayout->addWidget(m_mainWidget);

    m_mainLayout = new QVBoxLayout();
    m_mainWidget->setLayout(m_mainLayout);

    QLabel* positionSectionHeader = new QLabel(QString("Position of query point:"), this);
    positionSectionHeader->setStyleSheet("font-weight: bold");
    m_mainLayout->addWidget(positionSectionHeader);

    QGridLayout* positionLayout = new QGridLayout();

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

    QSpacerItem* buttonSpacer = new QSpacerItem(1, 20, QSizePolicy::Fixed, QSizePolicy::Fixed);
    m_mainLayout->addItem(buttonSpacer);

    m_liveUpdate = new QCheckBox(QString("Live Update"), this);
    m_mainLayout->addWidget(m_liveUpdate);

    QHBoxLayout* buttonLayout = new QHBoxLayout();

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

    connect(m_liveUpdate, SIGNAL(stateChanged(int)), this, SLOT(liveUpdateStateChanged(int)));

    connect(m_buttonApply, SIGNAL(pressed()), this, SIGNAL(applyPressed()));
    connect(m_buttonHide, SIGNAL(pressed()), this, SLOT(hidePress()));
}

void NearestNeighborWidget::deactivate()
{
    std::cerr << "Called deactivate in NearestNeighborWidget\n";
    m_mainWidget->setChecked(false);
}

void NearestNeighborWidget::activate()
{
    std::cerr << "Called activate in NearestNeighborWidget\n";
    m_mainWidget->setChecked(true);
}

void NearestNeighborWidget::resetValueRange(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax)
{
    m_xValueRange[0] = xMin;
    m_xValueRange[1] = xMax;

    m_yValueRange[0] = yMin;
    m_yValueRange[1] = yMax;

    m_zValueRange[0] = zMin;
    m_zValueRange[1] = zMax;


    updatePositionLabel();

    blockSliders(true);

    m_xValue->setValue(0);
    m_yValue->setValue(0);
    m_zValue->setValue(0);

    blockSliders(false);

    emit hidePressed();
    emit centerChanged(getXValue(), getYValue(), getZValue());
}

void NearestNeighborWidget::positionChanged(int)
{
    updatePositionLabel();
    emit centerChanged(getXValue(), getYValue(), getZValue());
}

void NearestNeighborWidget::liveUpdateStateChanged(int state)
{
    if (state == Qt::Unchecked)
    {
        emit liveUpdateChanged(false);
    }
    else
    {
        emit liveUpdateChanged(true);
    }
}

void NearestNeighborWidget::hidePress()
{
    m_liveUpdate->setChecked(false);
    emit hidePressed();
}

double NearestNeighborWidget::getXValue() const
{
    return getValue(m_xValue->value(), m_xValueRange);
}

double NearestNeighborWidget::getYValue() const
{
    return getValue(m_yValue->value(), m_yValueRange);
}

double NearestNeighborWidget::getZValue() const
{
    return getValue(m_zValue->value(), m_zValueRange);
}

double NearestNeighborWidget::getValue(int position, const double range[2]) const
{
    double factor = (static_cast<double>(position) - m_DEFAULT_SLIDER_MIN)
                    / (m_DEFAULT_SLIDER_MAX - m_DEFAULT_SLIDER_MIN);
    return (1.0 - factor)*range[0] + factor*range[1];
}

void NearestNeighborWidget::blockSliders(bool value)
{
    m_xValue->blockSignals(value);
    m_yValue->blockSignals(value);
    m_zValue->blockSignals(value);
}

void NearestNeighborWidget::updatePositionLabel()
{
    m_positionLabel->setText(QString("X: %1 \nY: %2 \nZ: %3").arg(QString::number(getXValue()),
                                                                  QString::number(getYValue()),
                                                                  QString::number(getZValue())));
}