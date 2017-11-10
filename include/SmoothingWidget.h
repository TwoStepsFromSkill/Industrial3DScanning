#ifndef SMOOTHING_WIDGET_H
#define SMOOTHING_WIDGET_H

#include "BaseTabWidget.h"

class QVBoxLayout;
class QPushButton;
class QDoubleSpinBox;
class QGroupBox;

class SmoothingWidget : public BaseTabWidget
{
    Q_OBJECT
public:
    SmoothingWidget(QWidget* parent = nullptr);
    virtual ~SmoothingWidget() {}

    virtual void deactivate() override;
    virtual void activate() override;

signals:
    void widgetEnabled(bool value);

    void applyPressed();

private:
    QVBoxLayout* m_dummyLayout;
    QGroupBox* m_mainWidget;

    QVBoxLayout* m_mainLayout;

    QDoubleSpinBox* m_radius;
    QPushButton* m_buttonApply;

    void blockSlider(bool value);
};

#endif // SMOOTHING_WIDGET_H
