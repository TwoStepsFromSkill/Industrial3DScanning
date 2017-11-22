#ifndef SMOOTHING_WIDGET_H
#define SMOOTHING_WIDGET_H

#include "BaseTabWidget.h"

class QVBoxLayout;
class QPushButton;
class QDoubleSpinBox;
class QGroupBox;
class QSlider;
class QRadioButton;
class QLabel;

class SmoothingWidget : public BaseTabWidget
{
    Q_OBJECT
public:
    SmoothingWidget(QWidget* parent = nullptr);
    virtual ~SmoothingWidget() {}

    virtual void deactivate() override;
    virtual void activate() override;

    void resetValueRange(double xMin, double xMax,
                       double yMin, double yMax,
                       double zMin, double zMax);
    void setNumberOfPoints(int number);

	void getRadius(double* radius);
    bool useGaussSmoothing() const;

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

    QGroupBox* m_smoothingMethodBox;
    QVBoxLayout* m_smoothingMethodLayout;
    QRadioButton* m_smoothMean;
    QRadioButton* m_smoothGauss;

    QPushButton* m_buttonApply;

    void blockSlider(bool value);

private slots:
    void radiusChanged(double val);
};

#endif // SMOOTHING_WIDGET_H
