include Make_linux.inc
#include Make_msys2.inc
#include Make_osx.inc

CXX = mpicxx	
CXXFLAGS = -std=c++17
ifeq ($(DEBUG),yes)
CXXFLAGS += -g -O0 -Wall -fbounds-check -pedantic -D_GLIBCXX_DEBUG
else
CXXFLAGS += -O3 -march=native -Wall
endif
CXXFLAGS += -fopenmp

CIBLE = ant_simu.exe
OBJETS = ant.o fractal_land.o renderer.o window.o ant_simu.o

LIB = -lpthread -lSDL2

default: help

all: $(CIBLE)

$(CIBLE): $(OBJETS)
	$(CXX) $(CXXFLAGS) $(OBJETS) -o $(CIBLE) $(LIB)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: all
	@mkdir -p results ../results
	@echo "Lancement de la simulation MPI avec 4 processus (GUI active)..."
	OMP_NUM_THREADS=2 mpirun -np 4 ./$(CIBLE) --gui 1 --iter 2000 --rep 1

benchmark: all
	@mkdir -p results ../results
	@echo "Exécution du benchmark MPI (4 processus, sans GUI)..."
	OMP_NUM_THREADS=2 mpirun -np 4 ./$(CIBLE) --gui 0 --iter 4000 --rep 5
	@cp results/Q4_timings_breakdown.csv ../results/Q4_mpi_timings.csv
	@echo "Résultats sauvegardés dans ../results/Q4_mpi_timings.csv"

run_scaling: all
	@mkdir -p results ../results
	@echo "Début de l'étude de scalabilité MPI (1, 2 et 4 processus)..."
	@for p in 1 2 4; do \
		echo "Exécution avec $$p processus MPI..."; \
		OMP_NUM_THREADS=2 mpirun -np $$p ./$(CIBLE) --gui 0 --iter 2000 --rep 3; \
		cp results/Q4_timings_breakdown.csv ../results/Q4_scaling_np$$p.csv; \
	done
	@echo "Étude terminée. Les fichiers Q4_scaling_npX.csv sont dans ../results/"

clean:
	@rm -f *.o *.exe *~

help:
	@echo "Cibles disponibles pour Q4 (MPI) :"
	@echo "  all          : Compilation de la version MPI"
	@echo "  run          : Exécution avec mpirun -np 4 (avec GUI)"
	@echo "  benchmark    : Mesure de performance avec 4 processus (sans GUI)"
	@echo "  run_scaling  : Test de performance de 1 à 5 processus"
	@echo "  clean        : Suppression des fichiers objets et de l'exécutable"
