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
MPI_BIND_FLAGS ?=
MPI_EXTRA_FLAGS ?=
P ?= 2
P_BENCH ?= 4
ITER ?= 4000
REP ?= 1
REP_BENCH ?= 5
ANTS ?= 5000
P_LIST ?= 1 2 4
LOCAL_RESULTS_DIR ?= results
RESULTS_DIR ?= ../results
CSV_NAME ?= Q5_timings_breakdown.csv

default: help

all: $(CIBLE)

$(CIBLE): $(OBJETS)
	$(CXX) $(CXXFLAGS) $(OBJETS) -o $(CIBLE) $(LIB)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: all
	@mkdir -p $(LOCAL_RESULTS_DIR) $(RESULTS_DIR)
	@echo "Running Q5 with P=$(P), ITER=$(ITER), REP=$(REP), ANTS=$(ANTS)"
	$(MPI_RUN) $(MPI_BIND_FLAGS) $(MPI_EXTRA_FLAGS) -np $(P) ./$(CIBLE) --gui 0 --iter $(ITER) --rep $(REP) --ants $(ANTS)
	@cp $(LOCAL_RESULTS_DIR)/$(CSV_NAME) $(RESULTS_DIR)/Q5_run_np$(P).csv

benchmark: all
	@mkdir -p $(LOCAL_RESULTS_DIR) $(RESULTS_DIR)
	@echo "Benchmark Q5 with P=$(P_BENCH), ITER=$(ITER), REP=$(REP_BENCH), ANTS=$(ANTS)"
	$(MPI_RUN) $(MPI_BIND_FLAGS) $(MPI_EXTRA_FLAGS) -np $(P_BENCH) ./$(CIBLE) --gui 0 --iter $(ITER) --rep $(REP_BENCH) --ants $(ANTS)
	@cp $(LOCAL_RESULTS_DIR)/$(CSV_NAME) $(RESULTS_DIR)/Q5_benchmark_p$(P_BENCH).csv
	@echo "Saved: $(RESULTS_DIR)/Q5_benchmark_p$(P_BENCH).csv"

run_scaling: all
	@mkdir -p $(LOCAL_RESULTS_DIR) $(RESULTS_DIR)
	@echo "Q5 scaling study for P_LIST=$(P_LIST)"
	@set -e; \
	for p in $(P_LIST); do \
		echo "Running np=$$p"; \
		$(MPI_RUN) $(MPI_BIND_FLAGS) $(MPI_EXTRA_FLAGS) -np $$p ./$(CIBLE) --gui 0 --iter $(ITER) --rep $(REP) --ants $(ANTS); \
		cp $(LOCAL_RESULTS_DIR)/$(CSV_NAME) $(RESULTS_DIR)/Q5_scaling_np$$p.csv; \
	done
	@echo "Scaling CSVs written to $(RESULTS_DIR)/Q5_scaling_npX.csv"

summarize_results:
	@mkdir -p $(RESULTS_DIR)
	@if [ ! -f "$(RESULTS_DIR)/Q5_scaling_np1.csv" ]; then \
		echo "Error: $(RESULTS_DIR)/Q5_scaling_np1.csv not found. Run make run_scaling first."; \
		exit 1; \
	fi
	@base_total=$$(awk -F, 'FNR==2 {print $$15}' "$(RESULTS_DIR)/Q5_scaling_np1.csv"); \
	out="$(RESULTS_DIR)/Q5_scaling_summary.csv"; \
	printf "P,T_compute_ms,T_comm_ms,T_total_ms,Speedup_vs_np1,Efficiency_vs_np1,Comm_pct,Ants_Min,Ants_Avg,Ants_Max,Imbalance\n" > "$$out"; \
	for f in $$(ls $(RESULTS_DIR)/Q5_scaling_np*.csv 2>/dev/null | sort -V); do \
		awk -F, -v base="$$base_total" 'FNR==2 {p=$$2; tc=$$10; tm=$$14; tt=$$15; sp=base/tt; eff=sp/p; cp=100*tm/tt; amin=$$16; aavg=$$17; amax=$$18; imb=$$20; printf "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.2f,%.2f,%.2f,%.2f,%.6f\n", p, tc, tm, tt, sp, eff, cp, amin, aavg, amax, imb }' "$$f" >> "$$out"; \
	done; \
	echo "Summary generated: $$out"; \
	cat "$$out"

clean:
	@rm -f *.o *.exe *~

help:
	@echo "Q5 targets:"
	@echo "  all                : Build Q5 MPI domain decomposition executable"
	@echo "  run                : One run with P=$(P)"
	@echo "  benchmark          : Benchmark run with P_BENCH=$(P_BENCH)"
	@echo "  run_scaling        : Scaling study for P_LIST='$(P_LIST)'"
	@echo "  summarize_results  : Summarize Q5_scaling_np*.csv"
	@echo "  clean              : Remove objects and executable"
	@echo "Variables: P, P_BENCH, ITER, REP, REP_BENCH, ANTS, P_LIST, MPI_BIND_FLAGS"
