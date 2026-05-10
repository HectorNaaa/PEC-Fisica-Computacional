% Ejercicio 2 — Patrones de Turing (Brusselator 2D)
% PEC Física Computacional
%
% Lee los datos generados por ejercicio1/main.c:
%   u_0.dat … u_9.dat  →  campo u(x,y) en 10 instantes
%   amplitud.dat        →  A(t) en cada paso temporal
%   time.dat            →  tiempos correspondientes
%
% Ejecutar desde la carpeta ejercicio2/

clear; clc; close all;

%% Parámetros
data_dir   = '../ejercicio1/';   % carpeta con los archivos .dat
L          = 100;
N_FRAMES   = 10;
DT         = 0.01;
FRAME_STEP = 2000;
frame_times = (0:N_FRAMES-1) * FRAME_STEP * DT;


%% Lectura de campos u_i.dat
fprintf('Cargando campos desde %s ...\n', data_dir);
U = cell(1, N_FRAMES);

for k = 0:N_FRAMES-1
    fname  = sprintf('u_%d.dat', k);
    U{k+1} = load(fullfile(data_dir, fname));
end


%% Escala de color global (fija para toda la animación)

umin = +Inf;
umax = -Inf;

for k = 1 : N_FRAMES
    umin = min(umin, min(U{k}(:)));
    umax = max(umax, max(U{k}(:)));
end

fprintf('Rango global: [%.4f, %.4f]\n\n', umin, umax);


%% Animación del patrón de Turing (Ejercicio 2.1)
% Frames iniciales: campo casi uniforme (ruido pequeño).
% Frames centrales: manchas emergentes con longitud de onda de Turing.
% Frames finales: patrón congelado por saturación no lineal.
% Causa: Dv >> Du → el inhibidor suprime al activador a distancia,
% pero el activador se autocataliza localmente → manchas.

fig_anim = figure('Name', 'Animación Patrón de Turing', ...
                  'Position', [100, 100, 600, 560]);

for k = 1 : N_FRAMES

    imagesc(U{k});
    colormap(parula);
    colorbar;
    caxis([umin, umax]);   % escala fija en todos los frames
    axis square;

    title(sprintf('Campo u(x,y) — Frame %d  (t = %.0f)', ...
                  k-1, frame_times(k)), 'FontSize', 13);
    xlabel('x'); ylabel('y');

    drawnow; pause(0.3);

end

% Guardar el último frame como PNG
last_frame_file = 'patron_turing_final.png';
saveas(fig_anim, last_frame_file);
fprintf('\nÚltimo frame guardado como: %s\n', last_frame_file);


%% Lectura de amplitud y tiempo (Ejercicio 2.2)
A = load(fullfile(data_dir, 'amplitud.dat'));
t = load(fullfile(data_dir, 'time.dat'));
fprintf('%d puntos temporales: t = [%.1f, %.1f]\n\n', length(t), t(1), t(end));


%% Representación semilogarítmica de A(t) (Ejercicio 2.2)
% A(t) = A0·exp(λ·t)  →  log(A) = log(A0) + λ·t  (lineal en semilogy).
% La pendiente de la zona recta es la tasa de crecimiento de Turing λ.

fig_amp = figure('Name', 'Evolución temporal de A(t)', ...
                 'Position', [750, 100, 700, 500]);

semilogy(t, A, 'b-', 'LineWidth', 1.2);
hold on;
xlabel('Tiempo t', 'FontSize', 12);
ylabel('A(t) = RMS de (u - u_0)', 'FontSize', 12);
title('Evolución de la amplitud de la perturbación', 'FontSize', 13);
grid on;


%% Ajuste lineal en el régimen exponencial
% Modificar t_ini/t_fin según la zona recta visible en semilogy.
t_ini = 20.0;
t_fin = 80.0;

idx_fit  = (t >= t_ini) & (t <= t_fin);
t_fit    = t(idx_fit);
logA_fit = log(A(idx_fit));

p      = polyfit(t_fit, logA_fit, 1);   % log(A) = p(1)*t + p(2)
lambda = p(1);   % tasa de crecimiento
logA0  = p(2);

A_ajustada = exp(polyval(p, t_fit));

% Añadir a la figura la recta de ajuste
semilogy(t_fit, A_ajustada, 'r--', 'LineWidth', 2.0);

legend({'A(t) simulación', sprintf('Ajuste: A_0·e^{\\lambda t},  \\lambda = %.4f', lambda)}, ...
       'Location', 'northwest', 'FontSize', 11);

% Marcar visualmente los límites del intervalo de ajuste
xline(t_ini, 'k:', 'LineWidth', 1.2);
xline(t_fin, 'k:', 'LineWidth', 1.2);

hold off;


%% Resultados
% λ > 0 → modo de Turing crece exponencialmente (inestabilidad confirmada).
% La longitud de onda dominante del patrón es L* = 2π/k* (ver Ejercicio 3).

fprintf('Ajuste en [%.1f, %.1f]:  lambda = %.6f,  A0 = %.4e\n', ...
        t_ini, t_fin, lambda, exp(logA0));
if lambda > 0
    fprintf('lambda > 0 -> inestabilidad de Turing confirmada.\n');
end
fprintf('\n');


%% Perfil espacial en la fila central

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


% Resultados esperados:
%   Frames 0-1: campo plano. Frames 3-6: manchas emergentes. Frames 7-9: saturado.
%   A(0) ≈ 0.01. Crecimiento exponencial t≈20-80. Saturación A≈0.28 (t>120).
%   lambda > 0 confirma inestabilidad de Turing para estos parámetros.
