#ifndef RANGE_QUERY_WIDGET_H
#define RANGE_QUERY_WIDGET_H

#include "BaseTabWidget.h"

class QVBoxLayout;
class QLabel;
class QPushButton;
class QSlider;
class QGroupBox;
class QCheckBox;

/**
 * @brief Widget that provides options to set up the range query.
 */
class RangeQueryWidget : public BaseTabWidget
{
    Q_OBJECT
public:
    /**
     * @brief Constructor.
     * @param parent Pointer to the parent element. (nullptr by default)
     */
    RangeQueryWidget(QWidget* parent = nullptr);
    /**
     * @brief Destructor.
     * @note Needs to be virtual since we are in a base class of an abstract
     *          class.
     */
    virtual ~RangeQueryWidget() {}

    /**
     * @brief Reset the range of possible coordiante values in x,y,z axis.
     *
     * This method has to called if the dataset changed. The widget needs
     * this information to adapt the range of the sliders to fit the dataset
     * extents.
     */
    void resetValueRange(double xMin, double xMax,
                       double yMin, double yMax,
                       double zMin, double zMax);

    /**
     * @brief Get the extents of the range query axis aligned box.
     * @param[in,out] minMax 6-element double array that is filled with the
     *                          extents of the range query box in the order
     *                          minX, maxX, minY, maxY, minZ, maxZ
     */
    void getBox(double* minMax);

    /**
     * @brief Deactivate the widget components.
     */
    virtual void deactivate() override;
    /**
     * @brief Activate the widget components.
     */
    virtual void activate() override;

signals:
    /**
     * @brief Signal to inform listeners that the widget was enable or disabled.
     */
    void widgetEnabled(bool value);

    /**
     * @brief Signal to inform listeners that the center of the box changed
     *          to x, y, z.
     */
    void centerChanged(double x, double y, double z);
    /**
     * @brief Signal to inform listeners that the box extents changed
     *          to dx, dy, dz.
     */
    void extendChanged(double dx, double dy, double dz);

    /**
     * @brief Signal to inform listeners that the live update option was toggled
     *          and has a new value.
     */
    void liveUpdateChanged(bool value);

    /**
     * @brief Signal to inform listeners that the apply button was pressed.
     */
    void applyPressed();
    /**
     * @brief Signal to inform listeners that the hide button was pressed.
     */
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
