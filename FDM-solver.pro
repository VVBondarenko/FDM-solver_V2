TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    poissontask.cpp \
    gnuplointerface.cpp \
    heattask.cpp

HEADERS += \
    poissontask.h \
    gnuplointerface.h \
    heattask.h

QMAKE_CXXFLAGS+= -fopenmp -fopenmp-simd
QMAKE_LFLAGS +=  -fopenmp -fopenmp-simd
