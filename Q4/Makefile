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

CIBLE = ant_simu.exe
OBJETS = ant.o fractal_land.o renderer.o window.o ant_simu.o

LIB = -lpthread -lSDL2

MPI_RUN ?= mpirun
MPI_BIND_FLAGS ?= --bind-to core --map-by core
MPI_EXTRA_FLAGS ?=
ITER ?= 4000
REP ?= 3
ANTS ?= 5000
P ?= 4
P_BENCH ?= 4
REP_BENCH ?= 5
P_LIST ?= 1 2 4
P_LIST_SMOKE ?= 8
LOCAL_RESULTS_DIR ?= results
RESULTS_DIR ?= ../results
CSV_NAME ?= Q4_timings_breakdown.csv

default: help

all: $(CIBLE)

$(CIBLE): $(OBJETS)
	$(CXX) $(CXXFLAGS) $(OBJETS) -o $(CIBLE) $(LIB)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: all
	@mkdir -p $(LOCAL_RESULTS_DIR) $(RESULTS_DIR)
	@echo "Execution MPI (P=$(P), binding coeur, sans GUI)..."
	$(MPI_RUN) $(MPI_BIND_FLAGS) $(MPI_EXTRA_FLAGS) -np $(P) ./$(CIBLE) --gui 0 --ants $(ANTS) --iter $(ITER) --rep 1
	@cp $(LOCAL_RESULTS_DIR)/$(CSV_NAME) $(RESULTS_DIR)/Q4_run_np$(P).csv

benchmark: all
	@mkdir -p $(LOCAL_RESULTS_DIR) $(RESULTS_DIR)
	@echo "Benchmark officiel MPI (P=$(P_BENCH), binding coeur, sans hwthreads)..."
	$(MPI_RUN) $(MPI_BIND_FLAGS) $(MPI_EXTRA_FLAGS) -np $(P_BENCH) ./$(CIBLE) --gui 0 --ants $(ANTS) --iter $(ITER) --rep $(REP_BENCH)
	@cp $(LOCAL_RESULTS_DIR)/$(CSV_NAME) $(RESULTS_DIR)/Q4_mpi_timings_p$(P_BENCH).csv
	@echo "Resultats sauvegardes dans $(RESULTS_DIR)/Q4_mpi_timings_p$(P_BENCH).csv"

run_scaling: all
	@mkdir -p $(LOCAL_RESULTS_DIR) $(RESULTS_DIR)
	@echo "Etude de scalabilite OFFICIELLE (P=$(P_LIST), binding coeur, sans hwthreads)."
	@echo "Si un P n'est pas disponible en coeurs reels, la commande echoue volontairement."
	@set -e; \
	for p in $(P_LIST); do \
		echo "Execution officielle avec $$p processus MPI..."; \
		$(MPI_RUN) $(MPI_BIND_FLAGS) $(MPI_EXTRA_FLAGS) -np $$p ./$(CIBLE) --gui 0 --ants $(ANTS) --iter $(ITER) --rep $(REP); \
		cp $(LOCAL_RESULTS_DIR)/$(CSV_NAME) $(RESULTS_DIR)/Q4_scaling_np$$p.csv; \
	done
	@echo "Etude officielle terminee. Fichiers: $(RESULTS_DIR)/Q4_scaling_npX.csv"

run_scaling_smoke: all
	@mkdir -p $(LOCAL_RESULTS_DIR) $(RESULTS_DIR)
	@echo "Smoke test MPI (non officiel) avec hwthreads pour P=$(P_LIST_SMOKE)."
	@set -e; \
	for p in $(P_LIST_SMOKE); do \
		echo "Smoke test avec $$p processus MPI (hwthreads)..."; \
		$(MPI_RUN) --use-hwthread-cpus -np $$p ./$(CIBLE) --gui 0 --ants $(ANTS) --iter $(ITER) --rep $(REP); \
		cp $(LOCAL_RESULTS_DIR)/$(CSV_NAME) $(RESULTS_DIR)/Q4_scaling_smoke_np$$p.csv; \
	done
	@echo "Smoke test termine. Fichiers: $(RESULTS_DIR)/Q4_scaling_smoke_npX.csv"

summarize_results:
	@mkdir -p $(RESULTS_DIR)
	@if [ ! -f "$(RESULTS_DIR)/Q4_scaling_np1.csv" ]; then \
		echo "Erreur: $(RESULTS_DIR)/Q4_scaling_np1.csv introuvable. Lancez d'abord make run_scaling."; \
		exit 1; \
	fi
	@base_total=$$(awk -F, 'FNR==2 {print $$13}' "$(RESULTS_DIR)/Q4_scaling_np1.csv"); \
	out="$(RESULTS_DIR)/Q4_scaling_summary.csv"; \
	printf "P,T_compute_ms,T_comm_ms,T_total_ms,T_comm_pher_ms,Speedup_vs_np1,Efficiency_vs_np1,Comm_pct\n" > "$$out"; \
	for f in $$(ls $(RESULTS_DIR)/Q4_scaling_np*.csv 2>/dev/null | sort -V); do \
		awk -F, -v base="$$base_total" 'FNR==2 {p=$$2; tc=$$9; tm=$$12; tt=$$13; tpher=$$11; sp=base/tt; eff=sp/p; cp=100*tm/tt; printf "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.2f\n", p, tc, tm, tt, tpher, sp, eff, cp }' "$$f" >> "$$out"; \
	done; \
	echo "Resume genere: $$out"; \
	cat "$$out"

clean:
	@rm -f *.o *.exe *~

help:
	@echo "Cibles disponibles pour Q4 (MPI) :"
	@echo "  all                 : Compilation de la version MPI pure"
	@echo "  run                 : Execution simple (P=$(P), binding coeur)"
	@echo "  benchmark           : Benchmark officiel (P=$(P_BENCH), sans hwthreads)"
	@echo "  run_scaling         : Scaling OFFICIEL pour P_LIST='$(P_LIST)'"
	@echo "  run_scaling_smoke   : Smoke test (non officiel) via --use-hwthread-cpus"
	@echo "  summarize_results   : Resume CSV des fichiers Q4_scaling_np*.csv"
	@echo "  clean               : Suppression des fichiers objets et de l'executable"
	@echo "Variables utiles: ITER, REP, ANTS, P, P_BENCH, REP_BENCH, P_LIST, P_LIST_SMOKE"
