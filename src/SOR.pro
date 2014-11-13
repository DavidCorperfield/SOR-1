TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

QMAKE_CXXFLAGS += -fopenmp -std=c++11
QMAKE_CXX = mpicxx

LIBS += -fopenmp -L/usr/local/openmpi/lib -lmpi_cxx -lstdc++ -ldl -lmpi -lblas

HEADERS += \
    boundaryconditions.hpp \
    grid.hpp \
    poissoneq.hpp \
    sor.hpp
