#ifndef THINNING_WIDGET_H
#define THINNING_WIDGET_H

#include "BaseTabWidget.h"

class QVBoxLayout;
class QPushButton;
class QGroupBox;
class QDoubleSpinBox;
class QSlider;
class QLabel;

class ThinningWidget : public BaseTabWidget
{
    Q_OBJECT
public:
    ThinningWidget(QWidget* parent = nullptr);
    virtual ~ThinningWidget() {}

    virtual void deactivate() override;
    virtual void activate() override;

    void resetValueRange(double xMin, double xMax,
                       double yMin, double yMax,
                       double zMin, double zMax);
    void setNumberOfPoints(int number);

    void getRadius(double* radius);

signals:
    void widgetEnabled(bool value);
    void tempPointChanged(int);

    void applyPressed();

private:
    QVBoxLayout* m_dummyLayout;
    QGroupBox* m_mainWidget;

    QVBoxLayout* m_mainLayout;

    QSlider* m_tmpPointSlider;

    QDoubleSpinBox* m_radius;
    QLabel* m_rangeLabel;

    QPushButton* m_buttonApply;

    void blockSlider(bool value);

private slots:
    void radiusChanged(double val);
};

#endif // THINNING_WIDGET_H
