include Make_linux.inc
#include Make_msys2.inc
#include Make_osx.inc

CXXFLAGS = -std=c++17
ifdef DEBUG
CXXFLAGS += -g -O0 -Wall -fbounds-check -pedantic -D_GLIBCXX_DEBUG
else
CXXFLAGS += -O3 -march=native -Wall
endif

CIBLE = ant_simu_q1.exe
OBJETS = ant.o fractal_land.o renderer.o window.o ant_simu.o

default: help

all: $(CIBLE)

$(CIBLE): $(OBJETS)
	$(CXX) $(CXXFLAGS) $(OBJETS) -o $(CIBLE) $(LIB)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: all
	@mkdir -p ../results
	@cd .. && ./Q1/$(CIBLE) --gui 0 --iter 3000 --rep 5

clean:
	@rm -f *.o *.exe *~

help:
	@echo "Cibles disponibles :"
	@echo "  all   : compilation de Q1"
	@echo "  run   : benchmark Q1 (CSV dans ../results)"
	@echo "  clean : nettoyage des artefacts"
	@echo "Variables utiles : DEBUG=yes"
