#include "SmoothingWidget.h"

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

SmoothingWidget::SmoothingWidget(QWidget* parent)
    : BaseTabWidget(parent)
{
    m_dummyLayout = new QVBoxLayout();
    m_dummyLayout->setContentsMargins(0, 0, 0, 0);
    m_mainWidget = new QGroupBox(this);
    m_mainWidget->setCheckable(true);
    m_mainWidget->setTitle(QString("Smoothing"));
    m_dummyLayout->addWidget(m_mainWidget);

    m_mainLayout = new QVBoxLayout();
    m_mainWidget->setLayout(m_mainLayout);

    m_radius = new QDoubleSpinBox(this);
    m_radius->setMinimum(0.0);
    m_radius->setSingleStep(0.01);
    m_radius->setMaximum(100.0);

    m_mainLayout->addWidget(m_radius);

    m_buttonApply = new QPushButton(QString("Apply"), this);
    m_mainLayout->addWidget(m_buttonApply);

    m_mainLayout->addStretch(1);

    setLayout(m_dummyLayout);

    connect(m_mainWidget, SIGNAL(toggled(bool)), this, SIGNAL(widgetEnabled(bool)));
    connect(m_buttonApply, SIGNAL(pressed()), this, SIGNAL(applyPressed()));

    deactivate();
}

void SmoothingWidget::deactivate()
{
    m_mainWidget->setChecked(false);
}

void SmoothingWidget::activate()
{
    m_mainWidget->setChecked(true);
}

void SmoothingWidget::blockSlider(bool value)
{
    m_radius->blockSignals(value);
}

void SmoothingWidget::getRadius(double* radius)
{
	*radius = m_radius->value();
}
