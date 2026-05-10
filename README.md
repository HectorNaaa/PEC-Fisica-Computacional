Disclaimer: Se ha usado como modelo de IA Claude Sonnet 4.6 desde visual studio code y git actions. Acerca del tema de la IA se ha usado en un 25-30% aproximadamente para limpieza, comentarios, revisión de los modelos numéricos implementados y ampliación y formateo del readme.me

# PEC — Física Computacional: Patrones de Turing

Simulación 2D del modelo Brusselator y análisis de la inestabilidad de Turing.

## Estructura

```
ejercicio1/   main.c                          — Simulación en C (Euler explícito)
ejercicio2/   ejercicio2.m                    — Visualización y ajuste λ (MATLAB)
ejercicio3/   ejercicio3_estabilidad_lineal.m — Estabilidad lineal (MATLAB)
```

## Uso

**1. Compilar y ejecutar la simulación (`ejercicio1/`)**

```bash
cd ejercicio1
gcc main.c -o simulacion -lm
./simulacion        # Windows: .\simulacion.exe
```

Genera: `u_0.dat … u_9.dat`, `amplitud.dat`, `time.dat`

**2. Análisis en MATLAB (`ejercicio2/`)**

```matlab
cd ejercicio2
run ejercicio2.m
```

Lee los `.dat` desde `../ejercicio1/`. Produce la animación del patrón y el ajuste de A(t) → λ.

**3. Estabilidad lineal (`ejercicio3/`)**

```matlab
cd ejercicio3
run ejercicio3_estabilidad_lineal.m
```

No depende de datos externos. Calcula el Jacobiano, autovalores y la relación de dispersión λ(k).

## Modelo

```
du/dt = Du·∇²u + a − u + u²v      Du = 1
dv/dt = Dv·∇²v + b − u²v          Dv = 10
```

Estado estacionario: `u0 = a+b = 1`,  `v0 = b/(a+b)² = 0.9`

La inestabilidad de Turing aparece porque `Dv >> Du`: el inhibidor v difunde mucho más rápido que el activador u, desestabilizando el estado homogéneo para modos espaciales con longitud de onda `L* = 2π/k*`. La λ* predicha en el ejercicio 3 debe coincidir con la pendiente medida empíricamente en el ejercicio 2.

Simulación numérica 2D del modelo Brusselator mediante Euler explícito.  
Genera patrones de Turing a partir de una perturbación gaussiana sobre el estado estacionario homogéneo.

---

## Conceptos físicos clave

### 1. ¿Qué significa físicamente el laplaciano discreto?

El laplaciano continuo `∇²u` mide la **curvatura local** del campo: si `u` en un punto es mayor que el promedio de sus vecinos, el laplaciano es negativo (el campo "cae" hacia los alrededores); si es menor, es positivo (el campo "sube"). En difusión, representa el **flujo neto de materia que entra o sale** de una celda por difusión.

En la red discreta con paso de red `dx = 1`, el laplaciano de los 5 puntos es:

```
Lap(u)[i,j] = u[i+1,j] + u[i-1,j] + u[i,j+1] + u[i,j-1] - 4·u[i,j]
```

Esto equivale a calcular cuánto difiere `u[i,j]` del **promedio de sus cuatro vecinos más próximos**, multiplicado por 4. Si el valor central es menor que ese promedio, la difusión lo aumenta; si es mayor, lo reduce. Este operador aproxima la segunda derivada espacial y converge al laplaciano continuo cuando `dx → 0`.

---

### 2. ¿Por qué se usan condiciones periódicas de contorno?

Las condiciones periódicas (o de contorno cíclicas) **conectan el borde derecho con el izquierdo y el borde superior con el inferior**, creando una topología toroidal. Esto tiene dos ventajas principales:

- **Elimina artefactos de frontera**: sin condiciones periódicas, el borde impone restricciones artificiales (Dirichlet, Neumann…) que distorsionan los patrones en el interior.
- **Simula un sistema efectivamente infinito**: todos los puntos de la red son equivalentes por simetría traslacional, lo que reproduce mejor el comportamiento de un sistema bulk.

En el código se implementan con aritmética modular:
```c
ip = (i + 1) % L;      // vecino derecho (vuelve al 0 cuando i = L-1)
im = (i - 1 + L) % L;  // vecino izquierdo (vuelve a L-1 cuando i = 0)
```

---

### 3. ¿Por qué el ruido inicial permite activar la inestabilidad de Turing?

El estado estacionario homogéneo `(u₀, v₀)` es un **punto fijo del sistema ODE** (sin difusión). Sin embargo, con difusión y parámetros adecuados (`Dv >> Du`), este estado puede ser inestable frente a **perturbaciones espaciales con cierta longitud de onda** — esto es la *inestabilidad de Turing* o *inestabilidad dirigida por difusión*.

El mecanismo es:
- El **activador** `u` difunde lentamente: si aumenta localmente, se autocataliza (`+u²v`).
- El **inhibidor** `v` difunde rápidamente: suprime el activador a larga distancia.
- El resultado es que ciertas longitudes de onda crecen exponencialmente.

Si partiéramos de `u[i,j] = u₀` exactamente (sin ruido), el sistema permanecería en el estado homogéneo **para siempre**, porque no hay ninguna perturbación que desencadene la inestabilidad. El **ruido gaussiano inicial** siembra todas las frecuencias espaciales simultáneamente con amplitud pequeña `σ`. La inestabilidad de Turing selecciona y amplifica las frecuencias dentro de la *banda inestable*, suprimiendo el resto, hasta que emergen los patrones visibles.

---

### 4. ¿Por qué se guarda A(t)?

La amplitud espacial RMS de la perturbación:

$$A(t) = \sqrt{\frac{1}{N} \sum_{i,j} \bigl(u_{i,j} - u_0\bigr)^2}$$

es el **observable más informativo** del proceso de formación del patrón:

| Régimen          | Comportamiento de A(t)                                |
|-----------------|------------------------------------------------------|
| Inicial          | A(t) ≈ σ = 0.01 (nivel del ruido)                   |
| Crecimiento lineal| A(t) ~ e^(λt), con λ = tasa de crecimiento de Turing|
| Saturación       | A(t) → constante (patrones no lineales maduros)      |

Guardando A(t) en cada paso se puede:
- **Verificar la inestabilidad**: la tasa de crecimiento exponencial debe coincidir con la predicción analítica del análisis de estabilidad lineal.
- **Identificar el tiempo de formación del patrón**.
- **Validar el método numérico**: si dt es demasiado grande, A(t) diverge (inestabilidad numérica).

---

## Estructura del código

```
Computational-Physics-PEC/
├── main.c          — Programa principal (Brusselator 2D)
├── README.md       — Este archivo
└── .gitignore
```

### Funciones implementadas

| Función                        | Descripción                                              |
|-------------------------------|----------------------------------------------------------|
| `gaussian_noise(sigma)`       | Genera ruido N(0,σ) mediante la transformada Box-Muller  |
| `initialize_fields(u0, v0)`   | Inicializa u y v con estado estacionario + ruido         |
| `compute_laplacian(field,i,j)`| Laplaciano discreto 5 puntos con contorno periódico      |
| `compute_amplitude(u0)`       | Calcula A(t) = RMS espacial de (u - u₀)                  |
| `save_field(frame)`           | Guarda u[i][j] en `u_<frame>.dat`                        |

---

## Compilación y ejecución (Windows, GCC)

Abre la **terminal de VS Code** (`` Ctrl+` ``) o PowerShell en la carpeta del proyecto:

```bash
# Compilar
gcc main.c -o simulacion -lm

# Ejecutar
./simulacion
```

> **Nota:** En Windows con MinGW puede que `./simulacion` no funcione directamente. Usa:
> ```bash
> .\simulacion.exe
> ```
> o simplemente:
> ```bash
> simulacion
> ```

---

## Archivos de salida

| Archivo              | Contenido                                              |
|---------------------|--------------------------------------------------------|
| `u_0.dat` … `u_9.dat`| Matriz L×L con el campo u(x,y) en 10 instantes        |
| `amplitud.dat`       | A(t) en cada paso temporal (20001 líneas)              |
| `time.dat`           | Tiempo t en cada paso (20001 líneas)                   |

Los snapshots se guardan cada 2000 pasos:

| Frame | Paso  | Tiempo  |
|-------|-------|---------|
| 0     | 0     | t = 0   |
| 1     | 2000  | t = 20  |
| …     | …     | …       |
| 9     | 18000 | t = 180 |

---

## Visualización con Python (opcional)

Puedes visualizar los patrones con:

```python
import numpy as np
import matplotlib.pyplot as plt

# Cargar un snapshot
u = np.loadtxt("u_9.dat")
plt.imshow(u, cmap="RdBu_r", origin="lower")
plt.colorbar(label="u(x,y)")
plt.title("Campo u — Patrón de Turing (t = 180)")
plt.tight_layout()
plt.savefig("turing_pattern.png", dpi=150)
plt.show()

# Graficar A(t)
t   = np.loadtxt("time.dat")
amp = np.loadtxt("amplitud.dat")
plt.semilogy(t, amp)
plt.xlabel("Tiempo t")
plt.ylabel("A(t)")
plt.title("Amplitud de la perturbación")
plt.tight_layout()
plt.savefig("amplitud.png", dpi=150)
plt.show()
```

---

## Parámetros de simulación

```
L     = 100      (red 100 × 100)
Du    = 1.0      (difusión activador)
Dv    = 10.0     (difusión inhibidor — Dv >> Du → inestabilidad)
a     = 0.1
b     = 0.9
dt    = 0.01
T     = 200.0
N     = 20000
sigma = 0.01

Estado estacionario:
  u0 = a + b = 1.0
  v0 = b / (a+b)^2 = 0.9
```

### Condición de estabilidad numérica (difusión)

Para que el método de Euler explícito sea estable para la difusión:

$$\Delta t \leq \frac{\Delta x^2}{2 D_{\max}} = \frac{1}{2 \times 10} = 0.05$$

Con `dt = 0.01 < 0.05`, el esquema es **estable**.
