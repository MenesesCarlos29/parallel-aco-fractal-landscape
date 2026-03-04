# parallel-aco-fractal-landscape

Guia rapida para compilar y ejecutar la simulacion ACO (Ant Colony Optimization) sobre paisaje fractal.

## 1. Requisitos

- Linux con `bash`
- Compilador `g++` con soporte C++17
- `make`
- SDL2 (para modo grafico `--gui 1`)

## 2. Compilacion

Desde la raiz del proyecto:

```bash
make -C src all
```

El ejecutable generado es:

- `src/ant_simu.exe`

## 3. Ejecucion rapida

### Modo grafico (por defecto)

```bash
./src/ant_simu.exe
```

### Modo benchmark (sin SDL)

```bash
./src/ant_simu.exe --gui 0 --iter 8000 --rep 5
```

Archivos generados:

- `results/Q0_baseline.csv`
- `results/machine_info.txt`

## 5. Parametros CLI del ejecutable

Comando base:

```bash
./src/ant_simu.exe [opciones]
```

| Parametro | Descripcion | Valor por defecto |
|---|---|---|
| `--grid-size N` | Tamano de la grilla `N x N` | `513` |
| `--ants N` | Numero de hormigas | `5000` |
| `--alpha X` | Coeficiente de ruido/propagacion de feromona (`0..1`) | `0.7` |
| `--beta X` | Coeficiente de evaporacion (`0..1`) | `0.999` |
| `--eps X` | Tasa de exploracion (`0..1`) | `0.8` |
| `--iter N` | Numero maximo de iteraciones por corrida | `1000` |
| `--rep N` | Numero de repeticiones (benchmark) | `1` |
| `--seed N` | Semilla inicial para reproducibilidad | `2026` |
| `--gui 0|1` | `1`: modo grafico, `0`: headless benchmark | `1` |
| `--help` | Muestra ayuda de parametros | N/A |

## 6. Notas importantes

- Si ejecutas con `--gui 1`, el programa fuerza internamente `--rep` a `1`.
- Para medir tiempo secuencial (`T_s`), usa siempre `--gui 0`.
- El CSV de Q0 incluye por repeticion: tiempo en ms, `Food_KPI` y `Checksum`, mas promedio y desviacion estandar al final.

## 7. Documentos del proyecto

- `PLAN_PROYECTO_Y_PROTOCOLO.md`: organizacion del trabajo, protocolo de benchmark y checklist.
- `ENUNCIADO_ACO_PARALELIZACION.md`: enunciado teorico completo del problema ACO y objetivos de paralelizacion.
- `Readme.pdf`: version PDF del enunciado original.
