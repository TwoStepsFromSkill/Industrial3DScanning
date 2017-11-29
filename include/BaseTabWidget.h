#ifndef BASE_TAB_WIDGET_H
#define BASE_TAB_WIDGET_H

#include <QWidget>

/**
 * @brief Base class widget for the widgets in the tab page.
 */
class BaseTabWidget : public QWidget
{
public:
    /**
     * @brief Constructor.
     * @param parent Pointer to parent widget.
     */
    BaseTabWidget(QWidget* parent) : QWidget(parent) {}

    /**
     * @brief Virtual destructor.
     */
    virtual ~BaseTabWidget() {}

    /**
     * @brief Deactivate widget. Will be overwritten be derived classes.
     */
    virtual void deactivate() = 0;

    /**
     * @brief Activate widget. Will be overwritten be derived classes.
     */
    virtual void activate() = 0;
};

#endif // BASE_TAB_WIDGET_H