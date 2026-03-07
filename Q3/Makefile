include Make_linux.inc
#include Make_msys2.inc
#include Make_osx.inc

CXXFLAGS = -std=c++17
ifeq ($(DEBUG),yes)
CXXFLAGS += -g -O0 -Wall -fbounds-check -pedantic -D_GLIBCXX_DEBUG
else
CXXFLAGS += -O3 -march=native -Wall
endif
CXXFLAGS += -fopenmp

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
	@./$(CIBLE) --gui 0 --iter 4000 --rep 5
	@if [ -f results/Q3_timings_breakdown.csv ]; then \
		cp results/Q3_timings_breakdown.csv ../results/Q3_timings_breakdown.csv; \
	else \
		echo "Erreur: fichier results/Q3_timings_breakdown.csv introuvable."; \
		exit 1; \
	fi
	@echo "Resultats enregistres dans ../results/Q3_timings_breakdown.csv"

run_q3: all
	@for t in 1 2 4 8; do \
		$(MAKE) --no-print-directory run_q3_thread THREADS=$$t; \
	done
	@echo "Resultats de scaling consolides dans ../results/Q3_openmp_scaling.csv"

run_q3_thread: all
	@mkdir -p results ../results
	@if [ -z "$(THREADS)" ]; then \
		echo "Erreur: variable THREADS manquante. Exemple: make run_q3_thread THREADS=2"; \
		exit 1; \
	fi
	@echo "Lancement Q3 avec OMP_NUM_THREADS=$(THREADS)"
	@OMP_NUM_THREADS=$(THREADS) ./$(CIBLE) --gui 0 --iter 4000 --rep 5
	@cp results/Q3_timings_breakdown.csv ../results/Q3_t$(THREADS)_timings_breakdown.csv
	@if [ "$(THREADS)" -eq 1 ]; then \
		cp results/Q3_timings_breakdown.csv ../results/Q3_threads1_reference.csv; \
	fi
	@mean_move=$$(awk -F, '$$1=="MEAN"{print $$5}' results/Q3_timings_breakdown.csv); \
	mean_evap=$$(awk -F, '$$1=="MEAN"{print $$6}' results/Q3_timings_breakdown.csv); \
	mean_update=$$(awk -F, '$$1=="MEAN"{print $$7}' results/Q3_timings_breakdown.csv); \
	mean_render=$$(awk -F, '$$1=="MEAN"{print $$8}' results/Q3_timings_breakdown.csv); \
	mean_total=$$(awk -F, '$$1=="MEAN"{print $$9}' results/Q3_timings_breakdown.csv); \
	std_total=$$(awk -F, '$$1=="STD_DEV"{print $$9}' results/Q3_timings_breakdown.csv); \
	food_rep1=$$(awk -F, 'NR==2{print $$3}' results/Q3_timings_breakdown.csv); \
	checksum_rep1=$$(awk -F, 'NR==2{print $$4}' results/Q3_timings_breakdown.csv); \
	if [ ! -f ../results/Q3_openmp_scaling.csv ]; then \
		echo "Threads,T_total_mean_ms,T_total_std_ms,T_move_ants_mean_ms,T_evap_mean_ms,T_pher_update_mean_ms,T_render_mean_ms,Food_KPI_rep1,Checksum_rep1" > ../results/Q3_openmp_scaling.csv; \
	fi; \
	awk -F, -v t="$(THREADS)" 'NR==1 || $$1!=t' ../results/Q3_openmp_scaling.csv > ../results/.Q3_openmp_scaling.tmp; \
	mv ../results/.Q3_openmp_scaling.tmp ../results/Q3_openmp_scaling.csv; \
	echo "$(THREADS),$$mean_total,$$std_total,$$mean_move,$$mean_evap,$$mean_update,$$mean_render,$$food_rep1,$$checksum_rep1" >> ../results/Q3_openmp_scaling.csv; \
	{ head -n 1 ../results/Q3_openmp_scaling.csv; tail -n +2 ../results/Q3_openmp_scaling.csv | sort -t, -k1,1n; } > ../results/.Q3_openmp_scaling.sorted; \
	mv ../results/.Q3_openmp_scaling.sorted ../results/Q3_openmp_scaling.csv
	@echo "Resultat mis a jour pour THREADS=$(THREADS) dans ../results/Q3_openmp_scaling.csv"

clean:
	@rm -f *.o *.exe *~

help:
	@echo "Cibles disponibles :"
	@echo "  all     : compilation de Q3 avec OpenMP"
	@echo "  run     : benchmark Q3 (CSV dans ../results/Q3_timings_breakdown.csv)"
	@echo "  run_q3  : benchmarks Q3 pour 1/2/4/8 threads (mise a jour du CSV global)"
	@echo "  run_q3_thread THREADS=N : benchmark Q3 pour un seul nombre de threads"
	@echo "  clean : nettoyage des artefacts"
	@echo "Variables utiles : DEBUG=yes"
