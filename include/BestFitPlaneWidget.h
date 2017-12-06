#ifndef BEST_FIT_PLANE_WIDGET_H
#define BEST_FIT_PLANE_WIDGET_H

#include "BaseTabWidget.h"

class QVBoxLayout;
class QPushButton;
class QGroupBox;

/**
 * @brief Widget that provides a button to start the best plane/line fit
 * algorithm.
 */
class BestFitPlaneWidget : public BaseTabWidget
{
    Q_OBJECT
public:
    /**
     * @brief Constructor.
     * @param parent Pointer to the parent element. (nullptr by default)
     */
    BestFitPlaneWidget(QWidget* parent = nullptr);

    /**
     * @brief Destructor.
     * @note Needs to be virtual since we are in a base class of an abstract
     *          class.
     */
    virtual ~BestFitPlaneWidget() {}

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
     * @brief Signal to inform listeners that the apply button was pressed.
     */
    void applyPressed();

private:
    QVBoxLayout* m_dummyLayout;
    QGroupBox* m_mainWidget;

    QVBoxLayout* m_mainLayout;

    QPushButton* m_buttonApply;
};

#endif // BEST_FIT_PLANE_WIDGET_H
