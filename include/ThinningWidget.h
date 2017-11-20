#ifndef THINNING_WIDGET_H
#define THINNING_WIDGET_H

#include "BaseTabWidget.h"

class QVBoxLayout;
class QPushButton;
class QGroupBox;
class QDoubleSpinBox;

class ThinningWidget : public BaseTabWidget
{
    Q_OBJECT
public:
    ThinningWidget(QWidget* parent = nullptr);
    virtual ~ThinningWidget() {}

    virtual void deactivate() override;
    virtual void activate() override;

    void getRadius(double* radius);

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

#endif // THINNING_WIDGET_H
