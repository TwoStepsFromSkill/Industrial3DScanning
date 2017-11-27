#include "ThinningWidget.h"
#include "Point3d.h"

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

ThinningWidget::ThinningWidget(QWidget* parent)
    : BaseTabWidget(parent)
{
    m_dummyLayout = new QVBoxLayout();
    m_dummyLayout->setContentsMargins(0, 0, 0, 0);
    m_mainWidget = new QGroupBox(this);
    m_mainWidget->setCheckable(true);
    m_mainWidget->setTitle(QString("Thinning"));
    m_dummyLayout->addWidget(m_mainWidget);

    m_mainLayout = new QVBoxLayout();
    m_mainWidget->setLayout(m_mainLayout);

    m_tmpPointSlider = new QSlider(this);
    m_tmpPointSlider->setMinimum(0);
    m_tmpPointSlider->setMaximum(1);
    m_tmpPointSlider->setTickInterval(1);
    m_tmpPointSlider->setOrientation(Qt::Horizontal);

    m_mainLayout->addWidget(m_tmpPointSlider);

    m_radius = new QDoubleSpinBox(this);
    m_radius->setMinimum(0.0);
    m_radius->setSingleStep(0.0001);
    m_radius->setDecimals(6);
    m_radius->setMaximum(100.0);

    m_mainLayout->addWidget(m_radius);

    m_rangeLabel = new QLabel(QString("..."), this);
    m_mainLayout->addWidget(m_rangeLabel);

    m_buttonApply = new QPushButton(QString("Apply"), this);
    m_mainLayout->addWidget(m_buttonApply);

    m_mainLayout->addStretch(1);

    setLayout(m_dummyLayout);

    connect(m_mainWidget, SIGNAL(toggled(bool)), this, SIGNAL(widgetEnabled(bool)));
    connect(m_buttonApply, SIGNAL(pressed()), this, SIGNAL(applyPressed()));

    connect(m_tmpPointSlider, SIGNAL(valueChanged(int)), this, SIGNAL(tempPointChanged(int)));
    connect(m_radius, SIGNAL(valueChanged(double)), this, SLOT(radiusChanged(double)));

    deactivate();
}

void ThinningWidget::resetValueRange(double xMin, double xMax, double yMin, double yMax,
                                      double zMin, double zMax)
{
    double diameter = distance3d(Point3d(xMin, yMin, zMin), Point3d(xMax, yMax, zMax));
    m_rangeLabel->setText(QString("Max diameter: %1").arg(QString::number(diameter)));
}

void ThinningWidget::setNumberOfPoints(int number)
{
    m_tmpPointSlider->setMinimum(0);
    m_tmpPointSlider->setMaximum(number - 1);
}

void ThinningWidget::deactivate()
{
    m_mainWidget->setChecked(false);
}

void ThinningWidget::activate()
{
    m_mainWidget->setChecked(true);
}

void ThinningWidget::blockSlider(bool value)
{
    m_radius->blockSignals(value);
}

void ThinningWidget::getRadius(double* radius)
{
	*radius = m_radius->value();
}

void ThinningWidget::radiusChanged(double)
{
    emit tempPointChanged(m_tmpPointSlider->value());
}
