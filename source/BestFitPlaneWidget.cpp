#include "BestFitPlaneWidget.h"

#include <QLabel>
#include <QString>
#include <QPushButton>
#include <QVBoxLayout>
#include <QGroupBox>

BestFitPlaneWidget::BestFitPlaneWidget(QWidget* parent)
    : BaseTabWidget(parent)
{
    m_dummyLayout = new QVBoxLayout();
    m_dummyLayout->setContentsMargins(0, 0, 0, 0);
    m_mainWidget = new QGroupBox(this);
    m_mainWidget->setCheckable(true);
    m_mainWidget->setTitle(QString("Best Fit Plane"));
    m_dummyLayout->addWidget(m_mainWidget);

    m_mainLayout = new QVBoxLayout();
    m_mainWidget->setLayout(m_mainLayout);

    m_buttonApply = new QPushButton(QString("Compute"), this);
    m_mainLayout->addWidget(m_buttonApply);

    m_mainLayout->addStretch(1);

    setLayout(m_dummyLayout);

    connect(m_mainWidget, SIGNAL(toggled(bool)), this, SIGNAL(widgetEnabled(bool)));
    connect(m_buttonApply, SIGNAL(pressed()), this, SIGNAL(applyPressed()));

    deactivate();
}

void BestFitPlaneWidget::deactivate()
{
    m_mainWidget->setChecked(false);
}

void BestFitPlaneWidget::activate()
{
    m_mainWidget->setChecked(true);
}
