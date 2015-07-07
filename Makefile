######################################################################
# Needed dependancies :                                              #
#                                                                    #
# - libfftw3-dev                                                     #
# - libsndfile1-dev                                                  #
# - libqt4-dev                                                       #
#                                                                    #
######################################################################

CXX		?= g++
MOC		?= moc
CXXFLAGS	+= -W -Wall -std=c++0x -pedantic
LDFLAGS		+= -lQtCore -lQtGui -lfftw3f -lsndfile
INCLUDE		+= -I /usr/include/qt4/QtCore -I /usr/include/qt4/
TARGET		= tadaridaD
TEST_SCRIPT = tests/test.sh

ifndef DEBUG
CXXFLAGS	+= -O2
else
CXXFLAGS	+= -g -ggdb -O0
endif

SRC		= main.cpp deteclaunch.cpp detec.cpp detectreatment.cpp
HEADERS		= deteclaunch.h detec.h
OBJ		= $(SRC:.cpp=.o)
MOC_OBJ		= $(HEADERS:.h=.moc.cpp)

all:$(TARGET)

$(TARGET): $(OBJ) $(MOC_OBJ)
	$(CXX) $^ -o $@ $(CXXFLAGS) $(INCLUDE) $(LDFLAGS)

%.moc.cpp: %.h
	$(MOC) $< -o $@

%.o:%.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS) $(INCLUDE)

clean:
	rm -f $(OBJ)
	rm -f $(TARGET)
	rm -f $(MOC_OBJ)
	rm -rf tests/waves/txt
	rm -rf log

test: all
	$(TEST_SCRIPT)
	$(TEST_SCRIPT) -p 6
	$(TEST_SCRIPT) -t 6
	$(TEST_SCRIPT) -t 3 -p 2
