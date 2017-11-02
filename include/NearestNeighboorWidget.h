#ifndef NEAREST_NEIGHBOOR_WIDGET_H
#define NEAREST_NEIGHBOOR_WIDGET_H

#include "BaseTabWidget.h"

class QVBoxLayout;
class QLabel;
class QPushButton;
class QSlider;
class QGroupBox;
class QCheckBox;

class NearestNeighboorWidget : public BaseTabWidget
{
    Q_OBJECT
public:
    NearestNeighboorWidget(QWidget* parent = nullptr);
    virtual ~NearestNeighboorWidget() {}

    void resetValueRange(double xMin, double xMax,
                       double yMin, double yMax,
                       double zMin, double zMax);

    virtual void deactivate() override;
    virtual void activate() override;

signals:
    void clickedEnable(bool value);

    void centerChanged(double x, double y, double z);

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

    QCheckBox* m_liveUpdate;

    QPushButton* m_buttonApply;
    QPushButton* m_buttonHide;

    double m_xValueRange[2];
    double m_yValueRange[2];
    double m_zValueRange[2];

    static constexpr int m_DEFAULT_SLIDER_MIN = -100;
    static constexpr int m_DEFAULT_SLIDER_MAX = 100;

    double getXValue() const;
    double getYValue() const;
    double getZValue() const;

    double getValue(int position, const double range[2]) const;

    void blockSliders(bool value);

    void updatePositionLabel();

private slots:
    void positionChanged(int);
    void liveUpdateStateChanged(int);

    void hidePress();
};

#endif // NEAREST_NEIGHBOOR_WIDGET_H
