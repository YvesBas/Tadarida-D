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

INCLUDEPATH += "C:/Program Files (x86)/Mega-Nerd/libsndfile/include" "Headers"

win32: LIBS += -L$$PWD/Libs/ -llibfftw3f-3 -llibsndfile-1

