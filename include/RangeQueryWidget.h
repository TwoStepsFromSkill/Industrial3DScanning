#ifndef RANGE_QUERY_WIDGET_H
#define RANGE_QUERY_WIDGET_H

#include "BaseTabWidget.h"

class QVBoxLayout;
class QLabel;
class QPushButton;
class QSlider;
class QGroupBox;
class QCheckBox;

class RangeQueryWidget : public BaseTabWidget
{
    Q_OBJECT
public:
    RangeQueryWidget(QWidget* parent = nullptr);
    virtual ~RangeQueryWidget() {}

    void resetValueRange(double xMin, double xMax,
                       double yMin, double yMax,
                       double zMin, double zMax);

    void getBox(double* minMax);

    virtual void deactivate() override;
    virtual void activate() override;

signals:
    void widgetEnabled(bool value);

    void centerChanged(double x, double y, double z);
    void extendChanged(double dx, double dy, double dz);

    void liveUpdateChanged(bool value);

    void applyPressed();
    void hidePressed();

private:
    QVBoxLayout* m_dummyLayout;
    QGroupBox* m_mainWidget;

    QVBoxLayout* m_mainLayout;

    QSlider* m_xValue;
    QSlider* m_yValue;
    QSlider* m_zValue;

    QLabel* m_positionLabel;

    QSlider* m_dxValue;
    QSlider* m_dyValue;
    QSlider* m_dzValue;

    QLabel* m_rangeLabel;

    QCheckBox* m_liveUpdate;

    QPushButton* m_buttonApply;
    QPushButton* m_buttonHide;

    double m_xValueRange[2];
    double m_yValueRange[2];
    double m_zValueRange[2];

    double m_dxValueRange[2];
    double m_dyValueRange[2];
    double m_dzValueRange[2];

    static constexpr int m_DEFAULT_SLIDER_MIN = -100;
    static constexpr int m_DEFAULT_SLIDER_MAX = 100;

    static constexpr int m_INITIAL_RANGE_VALUE = -75.0;

    double getXValue() const;
    double getYValue() const;
    double getZValue() const;

    double getXRange() const;
    double getYRange() const;
    double getZRange() const;

    double getValue(int position, const double range[2]) const;

    void blockSliders(bool value);

    void updatePositionLabel();
    void updateRangeLabel();

private slots:
    void positionChanged(int);
    void rangeChanged(int);
    void liveUpdateStateChanged(int);

    void hidePress();
};

#endif // RANGE_QUERY_WIDGET_H
