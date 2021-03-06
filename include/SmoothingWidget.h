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

/**
 * @brief Widget that provides options to set up the smoothing algorithm.
 */
class SmoothingWidget : public BaseTabWidget
{
    Q_OBJECT
public:
    /**
     * @brief Constructor.
     * @param parent Pointer to the parent element. (nullptr by default)
     */
    SmoothingWidget(QWidget* parent = nullptr);
    /**
     * @brief Destructor.
     * @note Needs to be virtual since we are in a base class of an abstract
     *          class.
     */
    virtual ~SmoothingWidget() {}

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
     * @brief Inform this widget about the number of points in the dataset.
     *
     * This information is required to provide a slider that can be used to
     * select every point in the dataset and displays its neighborhood
     * with the current settings.
     */
    void setNumberOfPoints(int number);

    /**
     * @brief Get the current value for the radius.
     * @param[in,out] Pointer to a double value that will be set to the radius.
     */
	void getRadius(double* radius);
    /**
     * @brief Check wether gaussian smoothing is enabled.
     */
    bool useGaussSmoothing() const;

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
     * @brief Signal to inform listeners that the temporary hint point changed.
     */
    void tempPointChanged(int);

    /**
     * @brief Signal to inform listeners that the apply button was pressed.
     */
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
