#ifndef NEAREST_NEIGHBOR_WIDGET_H
#define NEAREST_NEIGHBOR_WIDGET_H

#include "BaseTabWidget.h"

class QVBoxLayout;
class QLabel;
class QPushButton;
class QSlider;
class QGroupBox;
class QCheckBox;

/**
 * @brief Widget that provides options to set up the nearest neighbor search.
 */
class NearestNeighborWidget : public BaseTabWidget
{
    Q_OBJECT
public:
    /**
     * @brief Constructor.
     * @param parent Pointer to the parent element. (nullptr by default)
     */
    NearestNeighborWidget(QWidget* parent = nullptr);
    /**
     * @brief Destructor.
     * @note Needs to be virtual since we are in a base class of an abstract
     *          class.
     */
    virtual ~NearestNeighborWidget() {}

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
     * @brief Request location of query point for the nearest neighbor search.
     * @param[in,out] data Pointer to a 3-element double array. This will be
     *                      filled with the coordiantes of the query point.
     */
    void getQueryPoint(double* data);

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
     * @brief Signal to inform listeners that the query point location changed
     *          to x, y, z.
     */
    void positionChange(double x, double y, double z);

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

#endif // NEAREST_NEIGHBOR_WIDGET_H
