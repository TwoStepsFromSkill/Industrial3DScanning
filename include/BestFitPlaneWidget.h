#ifndef BEST_FIT_PLANE_WIDGET_H
#define BEST_FIT_PLANE_WIDGET_H

#include "BaseTabWidget.h"

class QVBoxLayout;
class QPushButton;
class QGroupBox;

class BestFitPlaneWidget : public BaseTabWidget
{
    Q_OBJECT
public:
    BestFitPlaneWidget(QWidget* parent = nullptr);
    virtual ~BestFitPlaneWidget() {}

    virtual void deactivate() override;
    virtual void activate() override;

signals:
    void widgetEnabled(bool value);

    void applyPressed();

private:
    QVBoxLayout* m_dummyLayout;
    QGroupBox* m_mainWidget;

    QVBoxLayout* m_mainLayout;

    QPushButton* m_buttonApply;
};

#endif // BEST_FIT_PLANE_WIDGET_H
