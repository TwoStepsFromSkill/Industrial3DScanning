#ifndef BASE_TAB_WIDGET_H
#define BASE_TAB_WIDGET_H

#include <QWidget>

class BaseTabWidget : public QWidget
{
public:
    BaseTabWidget(QWidget* parent) : QWidget(parent) {}
    virtual ~BaseTabWidget() {}

    virtual void deactivate() = 0;
    virtual void activate() = 0;
};

#endif // BASE_TAB_WIDGET_H