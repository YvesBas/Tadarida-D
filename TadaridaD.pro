#-------------------------------------------------
#
# Project created by QtCreator 2015-08-29T09:13:16
#
#-------------------------------------------------

QT       += core gui

QT       -= gui

TARGET = TadaridaD

TEMPLATE = app


SOURCES += main.cpp \
    deteclaunch.cpp \
    detec.cpp \
    detectreatment.cpp

HEADERS += \
    deteclaunch.h \
    detec.h \
    detectreatment.h

INCLUDEPATH += "C:/Program Files (x86)/Mega-Nerd/libsndfile/include" "Headers"

win32: LIBS += -L$$PWD/Libs/ -llibfftw3f-3 -llibsndfile-1
