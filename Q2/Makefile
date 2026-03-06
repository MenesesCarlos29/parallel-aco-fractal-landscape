include Make_linux.inc
#include Make_msys2.inc
#include Make_osx.inc

CXXFLAGS = -std=c++17
ifeq ($(DEBUG),yes)
CXXFLAGS += -g -O0 -Wall -fbounds-check -pedantic -D_GLIBCXX_DEBUG
else
CXXFLAGS += -O3 -march=native -Wall
endif

CIBLE = ant_simu.exe
OBJETS = ant.o fractal_land.o renderer.o window.o ant_simu.o

default: help

all: $(CIBLE)

$(CIBLE): $(OBJETS)
	$(CXX) $(CXXFLAGS) $(OBJETS) -o $(CIBLE) $(LIB)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: all
	@mkdir -p results ../results
	@./$(CIBLE) --gui 0 --iter 3000 --rep 5
	@if [ -f results/Q2_timings_breakdown.csv ]; then \
		cp results/Q2_timings_breakdown.csv ../results/Q2_timings_breakdown.csv; \
	else \
		echo "Erreur: fichier results/Q2_timings_breakdown.csv introuvable."; \
		exit 1; \
	fi
	@echo "Resultats enregistres dans ../results/Q2_timings_breakdown.csv"

clean:
	@rm -f *.o *.exe *~

help:
	@echo "Cibles disponibles :"
	@echo "  all   : compilation de Q2"
	@echo "  run   : benchmark Q2 (CSV dans ../results/Q2_timings_breakdown.csv)"
	@echo "  clean : nettoyage des artefacts"
	@echo "Variables utiles : DEBUG=yes"
