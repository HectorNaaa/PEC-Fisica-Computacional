% =========================================================================
% EJERCICIO 2 — Análisis de Patrones de Turing (Brusselator 2D)
% =========================================================================
% Física Computacional — PEC
%
% Este script analiza los datos generados por el programa en C:
%   u_0.dat … u_9.dat  →  campo u(x,y) en 10 instantes temporales
%   amplitud.dat        →  A(t) en cada paso temporal
%   time.dat            →  tiempos correspondientes
%
% Ejecutar desde la carpeta que contiene los archivos .dat
% =========================================================================

clear; clc; close all;

% =========================================================================
%% SECCIÓN 0 — PARÁMETROS DE LA SIMULACIÓN (para referencia y etiquetas)
% =========================================================================

L          = 100;       % Tamaño de la red (L x L)
N_FRAMES   = 10;        % Número de snapshots guardados
DT         = 0.01;      % Paso temporal del programa en C
FRAME_STEP = 2000;      % Cada cuántos pasos se guardó un frame

% Tiempo asociado a cada frame (t = frame * FRAME_STEP * DT)
frame_times = (0 : N_FRAMES-1) * FRAME_STEP * DT;
%   frame_times = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]


% =========================================================================
%% SECCIÓN 1 — LECTURA DE LOS CAMPOS u(x,y)
% =========================================================================
% Cargamos los 10 archivos en una estructura de celdas.
% Cada celda contiene una matriz L×L con el campo u en ese instante.

fprintf('=== Cargando archivos u_i.dat ===\n');
U = cell(1, N_FRAMES);       % Celda de matrices

for k = 0 : N_FRAMES - 1
    filename = sprintf('u_%d.dat', k);
    U{k+1}   = load(filename);   % load() detecta automáticamente matrices numéricas
    fprintf('  Leído: %s  (tamaño %dx%d)\n', filename, size(U{k+1},1), size(U{k+1},2));
end


% =========================================================================
%% SECCIÓN 2 — ESCALA DE COLOR GLOBAL
% =========================================================================
% Para comparar frames distintos necesitamos que la barra de color (caxis)
% sea la MISMA en todos ellos. Si usáramos el rango de cada frame por
% separado, los colores cambiarían de significado entre frames y la
% animación induciría a error.

umin = +Inf;
umax = -Inf;

for k = 1 : N_FRAMES
    umin = min(umin, min(U{k}(:)));
    umax = max(umax, max(U{k}(:)));
end

fprintf('\nRango global del campo u:\n');
fprintf('  umin = %.6f\n', umin);
fprintf('  umax = %.6f\n', umax);


% =========================================================================
%% SECCIÓN 3 — ANIMACIÓN DEL PATRÓN DE TURING (Ejercicio 2.1)
% =========================================================================
%
% Qué se observa físicamente:
% ─────────────────────────────────────────────────────────────────────────
% Al principio (frames 0-2) el campo es casi uniforme: el ruido inicial
% es tan pequeño (σ = 0.01) que apenas se distingue de la constante u₀=1.
%
% A partir del frame 3-5 (t ≈ 40-100) empiezan a aparecer MANCHAS o
% RAYAS con una longitud de onda característica λ_T, que es la del modo
% de Fourier más inestable según el análisis de estabilidad lineal.
%
% Este patrón persiste ("congela") en los últimos frames porque la
% no-linealidad satura el crecimiento: el sistema alcanza un nuevo
% estado estacionario INHOMOGÉNEO estable.
%
% Por qué aparecen las manchas:
% ─────────────────────────────────────────────────────────────────────────
% El activador u difunde lentamente (Du=1) pero se autocataliza (término u²v).
% El inhibidor v difunde rápidamente (Dv=10) y suprime a u a larga distancia.
% Esta separación de escalas (Dv >> Du) hace que el activador se concentre
% en islas rodeadas de zonas con inhibidor predominante → manchas de Turing.

fig_anim = figure('Name', 'Animación Patrón de Turing', ...
                  'Position', [100, 100, 600, 560]);

for k = 1 : N_FRAMES

    imagesc(U{k});                    % Muestra la matriz como imagen de color
    colormap(parula);                 % Mapa de colores continuo y perceptualmente uniforme
    colorbar;                         % Barra de escala de color a la derecha
    caxis([umin, umax]);              % Escala fija para todos los frames ← IMPORTANTE
    axis square;                      % Proporciones cuadradas (no distorsiona la red)

    title(sprintf('Campo u(x,y) — Frame %d  (t = %.0f)', ...
                  k-1, frame_times(k)), 'FontSize', 13);
    xlabel('x'); ylabel('y');

    drawnow;                          % Fuerza el redibujado inmediato en la ventana
    pause(0.3);                       % Espera 300 ms antes del siguiente frame

end

% Guardar el último frame como PNG
last_frame_file = 'patron_turing_final.png';
saveas(fig_anim, last_frame_file);
fprintf('\nÚltimo frame guardado como: %s\n', last_frame_file);


% =========================================================================
%% SECCIÓN 4 — LECTURA DE AMPLITUD Y TIEMPO (Ejercicio 2.2)
% =========================================================================

fprintf('\n=== Cargando amplitud.dat y time.dat ===\n');
A = load('amplitud.dat');   % Vector columna con A(t) en cada paso
t = load('time.dat');       % Vector columna con los tiempos

fprintf('  Número de puntos temporales: %d\n', length(t));
fprintf('  Tiempo inicial: %.3f  |  Tiempo final: %.3f\n', t(1), t(end));


% =========================================================================
%% SECCIÓN 5 — ESCALA SEMILOGARÍTMICA DE A(t)
% =========================================================================
%
% Por qué semilogy:
% ─────────────────────────────────────────────────────────────────────────
% En el régimen de inestabilidad lineal, la amplitud crece como:
%
%       A(t) = A₀ · exp(λ · t)
%
% Tomando logaritmo:
%
%       log(A(t)) = log(A₀) + λ · t     ← función LINEAL de t
%
% En escala semilogarítmica (eje y en log, eje x en lineal), una
% exponencial aparece como una RECTA PERFECTA. La pendiente de esa
% recta es exactamente λ, la tasa de crecimiento de Turing.
% Si la curva es cóncava (decrece) o convexa (se satura), estamos fuera
% del régimen lineal.

fig_amp = figure('Name', 'Evolución temporal de A(t)', ...
                 'Position', [750, 100, 700, 500]);

semilogy(t, A, 'b-', 'LineWidth', 1.2);
hold on;
xlabel('Tiempo t', 'FontSize', 12);
ylabel('A(t) = RMS de (u - u_0)', 'FontSize', 12);
title('Evolución de la amplitud de la perturbación', 'FontSize', 13);
grid on;


% =========================================================================
%% SECCIÓN 6 — AJUSTE LINEAL EN LA REGIÓN DE CRECIMIENTO EXPONENCIAL
% =========================================================================
%
% ─────────────────────────────────────────────────────────────────────────
% MODIFICAR AQUÍ el intervalo del ajuste si es necesario.
%   t_ini: tiempo donde empieza el crecimiento aproximadamente lineal
%   t_fin: tiempo donde la curva empieza a saturarse
% Elegir mirando la gráfica de semilogy: la zona recta es el régimen lineal.
% ─────────────────────────────────────────────────────────────────────────
t_ini = 20.0;    % Tiempo inicial del ajuste  (ajustar según la gráfica)
t_fin = 80.0;    % Tiempo final del ajuste    (ajustar según la gráfica)

% Seleccionar los índices dentro del intervalo [t_ini, t_fin]
idx_fit = (t >= t_ini) & (t <= t_fin);

t_fit    = t(idx_fit);          % Subvector de tiempos para el ajuste
A_fit    = A(idx_fit);          % Subvector de amplitudes para el ajuste
logA_fit = log(A_fit);          % log natural de A (lineariza la exponencial)

% Ajuste polinomial de grado 1:  logA = p(1)*t + p(2)
%   p(1) = λ  (pendiente = tasa de crecimiento)
%   p(2) = log(A₀)  (ordenada en el origen)
p = polyfit(t_fit, logA_fit, 1);

lambda = p(1);          % Tasa de crecimiento de Turing
logA0  = p(2);          % Logaritmo de la amplitud inicial ajustada

% Reconstruir la recta ajustada para representarla
A_ajustada = exp(polyval(p, t_fit));   % = exp(λ·t + log(A₀)) = A₀·exp(λ·t)

% Añadir a la figura la recta de ajuste
semilogy(t_fit, A_ajustada, 'r--', 'LineWidth', 2.0);

legend({'A(t) simulación', sprintf('Ajuste: A_0·e^{\\lambda t},  \\lambda = %.4f', lambda)}, ...
       'Location', 'northwest', 'FontSize', 11);

% Marcar visualmente los límites del intervalo de ajuste
xline(t_ini, 'k:', 'LineWidth', 1.2);
xline(t_fin, 'k:', 'LineWidth', 1.2);

hold off;


% =========================================================================
%% SECCIÓN 7 — RESULTADOS E INTERPRETACIÓN FÍSICA
% =========================================================================
%
% Significado físico de λ:
% ─────────────────────────────────────────────────────────────────────────
% λ es la TASA DE CRECIMIENTO del modo de Fourier más inestable.
%
%   λ > 0 → el modo crece exponencialmente → INESTABILIDAD DE TURING
%   λ < 0 → el modo decae → sistema estable (sin patrón)
%   λ = 0 → punto de bifurcación (umbral de Turing)
%
% En el análisis de estabilidad lineal alrededor del estado homogéneo,
% se busca la parte real máxima del espectro de la matriz jacobiana
% modificada por difusión:
%
%   M(k) = J - D·k²    donde J es el jacobiano cinético, D la matriz
%                        de difusión y k el número de onda.
%
% El modo k* con mayor Re[eigenvalue] de M(k*) domina el crecimiento:
% su tasa de crecimiento es exactamente λ. El patrón final tiene una
% longitud de onda espacial L* = 2π / k*.
%
% Por qué λ > 0 implica inestabilidad de Turing:
% ─────────────────────────────────────────────────────────────────────────
% El sistema sin difusión es estable (los eigenvalores del jacobiano J
% tienen parte real negativa). Cuando se añade difusión con Dv >> Du,
% el inhibidor se esparce rápidamente apagando el activador a distancia,
% pero localmente el activador tiene tiempo de autocatalizar antes de
% ser suprimido. Esto convierte eigenvalores negativos en positivos para
% un rango de k → inestabilidad que no existiría sin difusión.
%
% Relación entre λ y los modos de Fourier:
% ─────────────────────────────────────────────────────────────────────────
% El ruido inicial excita TODOS los modos de Fourier.
% Solo los modos con λ(k) > 0 crecen; los demás decaen.
% La amplitud total A(t) está dominada por el modo más inestable (λ_max),
% por eso observamos crecimiento exponencial con pendiente λ_max.
% La longitud de onda visible en los snapshots corresponde a k*.

fprintf('\n=================================================\n');
fprintf('  RESULTADOS — Ajuste lineal en [%.1f, %.1f]\n', t_ini, t_fin);
fprintf('=================================================\n');
fprintf('  Tasa de crecimiento:  lambda = %.6f\n', lambda);
fprintf('  Amplitud inicial:     A0     = %.6e\n', exp(logA0));
fprintf('\nInterpretación física:\n');
fprintf('  lambda > 0 (%.4f) confirma INESTABILIDAD DE TURING.\n', lambda);
fprintf('  El modo más inestable crece con e^(lambda·t).\n');
fprintf('  La longitud de onda dominante del patrón final es\n');
fprintf('  la del modo k* que maximiza la tasa de crecimiento.\n');
fprintf('=================================================\n\n');


% =========================================================================
%% SECCIÓN 8 — FIGURA ADICIONAL: PERFIL DE u A LO LARGO DE UNA FILA
% =========================================================================
% Útil para estimar visualmente la longitud de onda espacial del patrón.

fig_perfil = figure('Name', 'Perfil espacial u(x) — último frame', ...
                    'Position', [100, 700, 700, 300]);

fila = round(L / 2);                        % Fila central de la red
plot(1:L, U{N_FRAMES}(fila, :), 'b-', 'LineWidth', 1.5);
xlabel('x (índice de celda)', 'FontSize', 12);
ylabel('u(x, L/2)', 'FontSize', 12);
title(sprintf('Perfil de u a lo largo de y = %d  (frame %d, t = %.0f)', ...
              fila, N_FRAMES-1, frame_times(N_FRAMES)), 'FontSize', 12);
grid on;

% Guardar figura de amplitud
saveas(fig_amp,    'amplitud_turing.png');
saveas(fig_perfil, 'perfil_turing.png');
fprintf('Figuras adicionales guardadas:\n');
fprintf('  amplitud_turing.png\n');
fprintf('  perfil_turing.png\n');


% =========================================================================
%% SECCIÓN 9 — RESULTADOS ESPERADOS (referencia para corrección)
% =========================================================================
%
% Si la simulación está bien implementada, se deberían observar:
%
%  1. CAMPO ESPACIAL (Ej. 2.1):
%     - Frames 0-1: campo casi uniforme ~ u₀ = 1.0, sin estructura visible.
%     - Frames 3-5: aparición de manchas irregulares con una escala espacial
%       característica (longitud de onda de Turing).
%     - Frames 7-9: patrón "congelado" estacionario con manchas bien definidas
%       de activador (azul claro) sobre fondo de inhibidor (azul oscuro).
%       El patrón NO debe cambiar de aspecto entre frames 8 y 9 → saturación.
%
%  2. AMPLITUD A(t) (Ej. 2.2):
%     - A(0) ≈ σ = 0.01 (nivel del ruido inicial).
%     - Crecimiento exponencial claro en escala semilogarítmica entre
%       t ≈ 20 y t ≈ 80, con pendiente λ > 0 (λ ≈ 0.03 – 0.10 típico).
%     - Saturación: A(t) → constante ≈ 0.28 para t > 120.
%     - El ajuste lineal en el régimen exponencial debe dar R² ≈ 1.
%
%  3. INTERPRETACIÓN:
%     - λ > 0 confirma la inestabilidad de Turing para estos parámetros.
%     - La longitud de onda dominante puede estimarse del perfil espacial
%       como la distancia media entre manchas consecutivas.
%     - Si λ fuera negativo, no habría formación de patrón.
%
% =========================================================================
