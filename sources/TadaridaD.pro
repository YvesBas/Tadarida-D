#-------------------------------------------------
#
# Project created by QtCreator 2014-12-16T17:17:09
#
#-------------------------------------------------

win32 {
DEFINES += BUILDTIME=\\\"$$system('echo %time%')\\\"
DEFINES += BUILDDATE=\\\"$$system('echo %date%')\\\"
} else {
DEFINES += BUILDTIME=\\\"$$system(date '+%H:%M.%s')\\\"
DEFINES += BUILDDATE=\\\"$$system(date '+%d/%m/%y')\\\"
}

QT       += core


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

INCLUDEPATH += "Headers"

LIBS += -L$$PWD/Libs/ -lfftw3f -lsndfile


