VPATH = ../src
CC = $(CXX)

ALLFLAGS = -std=gnu++17 -Wall -O3 -fopenmp -march=native
#ALLFLAGS = -std=gnu++17 -Wall -O3 -march=native

CXXFLAGS = $(PROGVAR) $(ALLFLAGS) $(shell sdl2-config --cflags) -I/usr/include/eigen3
LDFLAGS = $(ALLFLAGS)
LDLIBS = -lstdc++fs $(shell sdl2-config --libs)

objects = agents.o agent.o agentinterface.o agentmanager.o calibrator.o colorfinder.o diffuser.o diffuserfdm.o fluidvel.o managerinterface.o picture.o screen.o nucleus.o world.o

agents : $(objects)

depend :
	g++ -MM ../src/*.cpp > dependfile

.PHONY : clean
clean :
	rm agents $(objects)

include dependfile
