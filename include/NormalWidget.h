#ifndef NORMAL_WIDGET_H
#define NORMAL_WIDGET_H

#include "BaseTabWidget.h"

class QVBoxLayout;
class QPushButton;
class QGroupBox;
class QDoubleSpinBox;
class QSlider;
class QLabel;

class NormalWidget : public BaseTabWidget
{
    Q_OBJECT
public:
    /**
     * @brief Constructor.
     * @param parent Pointer to the parent element. (nullptr by default)
     */
    NormalWidget(QWidget* parent = nullptr);
    /**
     * @brief Destructor.
     * @note Needs to be virtual since we are in a base class of an abstract
     *          class.
     */
    virtual ~NormalWidget() {}

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
     * @brief Get the current value for the radius.
     * @param[in,out] Pointer to a double value that will be set to the radius.
     */
	void getRadius(double* radius);

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

    QDoubleSpinBox* m_radius;
    QLabel* m_rangeLabel;

    QPushButton* m_buttonApply;
};

#endif // NORMAL_WIDGET_H
