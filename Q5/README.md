# Q5 - MPI enfoque 2 (Domain Decomposition por filas)

## 1) Resumen de implementacion

`Q5` implementa la simulacion ACO en **MPI puro** con **descomposicion 1D por filas** (row-wise):

- Cada rank es dueno de un bloque contiguo de filas globales `[row_start, row_end)`.
- Cada rank mantiene solo su subdominio local de:
  - terreno (`fractal_land`),
  - feromonas (`pheronome`),
  - mas dos filas halo (ghost rows): superior e inferior.
- Las hormigas pertenecen al rank que posee su coordenada global `y`.
- Si una hormiga cruza el borde del subdominio, migra al rank vecino (arriba/abajo).
- La coherencia de feromonas entre subdominios se mantiene con **halo exchange entre vecinos**.

No se usa entorno replicado ni reduccion global del mapa de feromonas.

## 2) Estrategia de dominio y comunicacion

### 2.1 Particion del dominio

- Dominio global: `N x N`.
- Reparto por filas entre `P` ranks:
  - `row_start(rank)`,
  - `local_height(rank)`,
  - `row_end(rank) = row_start + local_height`.

### 2.2 Ghost cells / halo rows

Cada rank almacena `local_height + 2` filas:

- fila local `1..local_height` (interior),
- fila `0` (halo superior),
- fila `local_height+1` (halo inferior).

Las ghost rows se usan para que una hormiga en el borde pueda consultar feromonas vecinas entre ranks sin reconstruir mapa global.

### 2.3 Halo exchange

En cada iteracion, despues de actualizar feromonas:

1. envio/recepcion con vecino superior,
2. envio/recepcion con vecino inferior.

Se hace con `MPI_Sendrecv` entre vecinos directos.
No hay `MPI_Allreduce` del mapa completo.

### 2.4 Migracion de hormigas

Cuando una hormiga sale del rango de filas local:

- se empaqueta y envia al vecino correspondiente,
- se recibe migracion entrante,
- se conserva estado completo (`x`, `y`, `is_loaded`, `seed`, `consumed_time`, `id`).

La conservacion del total de hormigas se valida con `Total_Ants_Global`.

## 3) Metricas y validacion

### 3.1 Metricas temporales (MPI_Wtime)

Se reportan (en ms):

- `T_compute_move_ms`
- `T_compute_evap_ms`
- `T_compute_update_ms`
- `T_compute_ms`
- `T_comm_halo_ms`
- `T_comm_migration_ms`
- `T_comm_food_ms`
- `T_comm_ms`
- `T_total_ms`

Los tiempos globales reportados en rank 0 se agregan con `MPI_Reduce(..., MPI_MAX, ...)`.

### 3.2 Metricas funcionales

- `Food_KPI` global.
- `Checksum` global determinista en rank 0 (gather de estado global).
- `Total_Ants_Global` (conservacion de hormigas).
- Carga por rank:
  - `Ants_Min`
  - `Ants_Avg`
  - `Ants_Max`
  - `Ants_Stddev`
  - `Imbalance_Ratio = max / avg`.

### 3.3 Validacion contra Q2

- Para `np=1`, se compara contra `Q2/results/Q2_timings_breakdown.csv`.
- La verificacion se reporta explicitamente como:
  - `MATCH_VALUE_ONLY`,
  - `MISMATCH`,
  - o `A_VERIFIER` con razon.

## 4) Compilacion y ejecucion

Desde raiz del repo:

```bash
cd Q5
make clean && make all
```

Ejecucion simple:

```bash
make run P=2 ITER=200 REP=1 ANTS=5000
```

Benchmark:

```bash
make benchmark P_BENCH=4 ITER=2000 REP_BENCH=5 ANTS=5000
```

Scaling:

```bash
make run_scaling ITER=200 REP=1 ANTS=5000 P_LIST="1 2 4"
make summarize_results
```

Ejecucion manual:

```bash
mpirun -np 4 ./Q5/ant_simu.exe --gui 0 --iter 200 --rep 1 --ants 5000
```

`rank 0` es el unico que escribe CSV.

## 5) Resultados observados (smoke de referencia)

Configuracion usada: `ITER=200`, `REP=1`, `ANTS=5000`, `P in {1,2,4}`.

| P | T_compute (ms) | T_comm (ms) | T_total (ms) | Speedup vs P=1 | Efficiency | Comm % |
|---|----------------|-------------|--------------|----------------|------------|--------|
| 1 | 623.904 | 1.871 | 624.664 | 1.000 | 1.000 | 0.30% |
| 2 | 339.261 | 45.985 | 348.466 | 1.793 | 0.896 | 13.20% |
| 4 | 210.221 | 78.610 | 227.483 | 2.746 | 0.686 | 34.56% |

Observacion clave: sube `T_comm` con `P`, aunque `T_compute` baja.

## 6) Analisis tecnico de escalabilidad

### Por que puede degradarse al subir `np`

- Halo exchange ocurre en cada iteracion y crece el costo de sincronizacion.
- La migracion de hormigas agrega trafico adicional entre vecinos.
- Aunque se reparte el movimiento de hormigas, crece la fraccion de tiempo de comunicacion.
- El nido/comida puede concentrar hormigas en ciertas franjas de filas, aumentando desbalance.
- Con desbalance, el tiempo real lo marca el rank mas lento (reduccion por `MPI_MAX`).

Esto es coherente con el enfoque de domain decomposition y no implica por si solo error de medicion.

## 7) Limitaciones actuales

- La referencia `Q2` para `np=1` actualmente reporta `MISMATCH` en checksum (se reporta de forma honesta).
- El balanceo es estatico por filas; no hay rebalanceo dinamico.
- Con cargas pequenas, la comunicacion puede dominar antes de observar mejor escalado.

## 8) Mitigaciones dentro del alcance de Q5

- Ejecutar con tamanos/problemas mayores para amortizar comunicacion.
- Afinar `P` segun hardware real para evitar sobrerreparto.
- Mantener medicion separada `compute` vs `comm` y monitorear `Imbalance_Ratio`.

## 9) Trabajo futuro (fuera de alcance de Q5)

- Rebalanceo dinamico de dominio.
- Migracion no bloqueante/pipeline mas avanzado.
- Esquemas de comunicacion mas selectivos para reducir trafico.
