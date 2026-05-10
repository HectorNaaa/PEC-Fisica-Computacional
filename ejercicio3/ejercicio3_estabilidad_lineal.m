% =========================================================================
% PEC 2026 - Ejercicio 3
% Análisis de estabilidad lineal del modelo reacción-difusión (Brusselator)
% Este script no requiere datos externos.
% Autor: HectorNaaa
% =========================================================================
%
% Modelo:
%   du/dt = Du·∇²u + f(u,v)    f(u,v) = a - u + u²v
%   dv/dt = Dv·∇²v + g(u,v)    g(u,v) = b - u²v
%
% El análisis de estabilidad lineal responde a la pregunta:
%   "¿Sobreviven las perturbaciones pequeñas sobre el estado estacionario
%    homogéneo (u0, v0)?"
%
% Sin difusión → Estabilidad local (ODE)
% Con difusión → Posible inestabilidad de Turing (PDE)
% =========================================================================

clear; clc; close all;

%% Parámetros
a  = 0.1;
b  = 0.9;
Du = 1.0;
Dv = 10.0;   % Dv >> Du: condición clave para la inestabilidad de Turing

% Estado estacionario: f(u0,v0)=0, g(u0,v0)=0  →  u0=a+b, v0=b/u0²
u0 = a + b;
v0 = b / u0^2;

fprintf('Parámetros: a=%.2f  b=%.2f  Du=%.1f  Dv=%.1f\n', a, b, Du, Dv);
fprintf('Estado estacionario: u0=%.4f  v0=%.4f\n\n', u0, v0);


%% SECCIÓN 1 — JACOBIANO DEL SISTEMA SIN DIFUSIÓN
%
% Linealizando alrededor de (u0, v0):
%   f = a - u + u²v  →  ∂f/∂u = -1+2uv,  ∂f/∂v = u²
%   g = b - u²v      →  ∂g/∂u = -2uv,    ∂g/∂v = -u²

J11 = -1 + 2*u0*v0;   % ∂f/∂u
J12 =  u0^2;           % ∂f/∂v
J21 = -2*u0*v0;        % ∂g/∂u
J22 = -u0^2;           % ∂g/∂v

J = [J11, J12;
     J21, J22];

fprintf('--- SECCIÓN 1: Jacobiano sin difusión ---\n\n');
fprintf('J = [ %+.4f   %+.4f ]\n', J(1,1), J(1,2));
fprintf('    [ %+.4f   %+.4f ]\n\n', J(2,1), J(2,2));

% Propiedades de J que determinan la estabilidad:
traza_J    = trace(J);   % tr(J) = J11 + J22 < 0 es condición necesaria de estabilidad
det_J      = det(J);     % det(J) > 0 es condición necesaria de estabilidad
fprintf('tr(J)  = %.6f  (debe ser < 0 para estabilidad)\n', traza_J);
fprintf('det(J) = %.6f  (debe ser > 0 para estabilidad)\n\n', det_J);


%% SECCIÓN 2 — AUTOVALORES SIN DIFUSIÓN
%
% Perturbación homogénea ~ exp(λ·t).
% Re(λ) < 0 para ambos → estado estable. Complejos con Re<0 → foco estable.

lambda_sin_difusion = eig(J);   % MATLAB devuelve los 2 autovalores del sistema 2x2

fprintf('--- SECCIÓN 2: Autovalores sin difusión ---\n\n');
for i = 1:length(lambda_sin_difusion)
    lam = lambda_sin_difusion(i);
    fprintf('  λ_%d = %+.6f %+.6fi    Re(λ_%d) = %+.6f\n', ...
            i, real(lam), imag(lam), i, real(lam));
end
fprintf('\n');

partes_reales = real(lambda_sin_difusion);

if all(partes_reales < 0)
    fprintf('CONCLUSIÓN: Re(λ) < 0 para todos los autovalores.\n');
    fprintf('→ El estado estacionario es ESTABLE sin difusión.\n');
    if ~isreal(lambda_sin_difusion)
        fprintf('→ Los autovalores son complejos: el punto fijo es un FOCO ESTABLE.\n');
        fprintf('  (Las perturbaciones homogéneas oscilan mientras decaen.)\n');
    else
        fprintf('→ Los autovalores son reales: el punto fijo es un NODO ESTABLE.\n');
    end
else
    fprintf('ADVERTENCIA: Re(λ) > 0 para algún autovalor.\n');
    fprintf('→ El estado estacionario es INESTABLE incluso sin difusión.\n');
end

% Inestabilidad de Turing: el sistema es estable sin difusión pero
% se desestabiliza al añadir difusión con Dv >> Du.
fprintf('\n');


%% SECCIÓN 3 — RELACIÓN DE DISPERSIÓN CON DIFUSIÓN
%
% Con difusión, una perturbación ~ exp(λ·t)·cos(k·x) evoluciona según:
%   M(k) = J - k²·diag(Du, Dv)
% Si Re(λ_dom(k)) > 0 para algún k, hay inestabilidad de Turing.

k_values = linspace(0, 2, 1000);    % Vector de números de onda (rad/píxel)
N_k = length(k_values);
lambda_values = zeros(1, N_k);      % Almacenará λ_max(k) = max(Re(eig(M(k))))

for idx = 1:N_k
    k  = k_values(idx);
    k2 = k^2;

    % Matriz de estabilidad modificada por difusión:
    M = [J11 - Du*k2,   J12;
         J21,            J22 - Dv*k2];

    % Autovalores de M(k) — el dominante es el de mayor parte real
    eigs_k = eig(M);
    lambda_values(idx) = max(real(eigs_k));
end


%% SECCIÓN 4 — MODO DOMINANTE Y ESCALA ESPACIAL
%
% k* = argmax λ(k): número de onda con mayor tasa de crecimiento.
% L* = 2π/k*: longitud de onda dominante del patrón emergente.

[lambda_max, idx_max] = max(lambda_values);
k_max = k_values(idx_max);

if k_max > 0
    longitud_onda = 2*pi / k_max;
else
    longitud_onda = Inf;
end

fprintf('--- SECCIÓN 3+4: Relación de dispersión con difusión ---\n\n');
fprintf('  k*        = %.6f  rad/píxel\n', k_max);
fprintf('  λ_max(k*) = %.6f\n', lambda_max);
fprintf('  L* = 2π/k* = %.4f  píxeles (longitud de onda dominante)\n', longitud_onda);
fprintf('\n');

if lambda_max > 0
    fprintf('CONCLUSIÓN: λ_max = %.6f > 0\n', lambda_max);
    fprintf('→ Existe una banda de números de onda k con λ(k) > 0.\n');
    fprintf('→ INESTABILIDAD DE TURING confirmada.\n');
    fprintf('→ El patrón emergente tendrá una escala espacial ≈ %.2f píxeles.\n', longitud_onda);
    fprintf('→ Compare con la longitud de onda visual en los snapshots del Ej.1.\n');
else
    fprintf('CONCLUSIÓN: λ_max = %.6f ≤ 0\n', lambda_max);
    fprintf('→ No hay inestabilidad de Turing para estos parámetros.\n');
end
fprintf('\n');


%% SECCIÓN 5 — FIGURA: RELACIÓN DE DISPERSIÓN
% Banda de k con λ(k) > 0: modos inestables. k* da el modo dominante.

fig = figure('Name', 'Relación de dispersión — Brusselator', ...
             'Position', [100, 100, 800, 500]);

% Curva principal
plot(k_values, lambda_values, 'b-', 'LineWidth', 2.0, ...
     'DisplayName', '\lambda_{max}(k)');
hold on;

% Línea de referencia λ = 0 (umbral de estabilidad)
yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', '\lambda = 0  (umbral)');

% Marcar el máximo k*
plot(k_max, lambda_max, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', ...
     'DisplayName', sprintf('k^* = %.4f,  \\lambda_{max} = %.4f', k_max, lambda_max));

% Anotación del máximo
annotation_text = sprintf('  k^* = %.3f\n  \\lambda = %.4f\n  L^* = %.2f', ...
                           k_max, lambda_max, longitud_onda);
text(k_max + 0.04, lambda_max, annotation_text, ...
     'FontSize', 10, 'Color', 'r');

hold off;

xlabel('k  (número de onda, rad/píxel)', 'FontSize', 12);
ylabel('Re(\lambda_{dominante})',          'FontSize', 12);
title('Relación de dispersión del sistema reacción-difusión (Brusselator)', ...
      'FontSize', 13);
legend('Location', 'northeast', 'FontSize', 11);
grid on;

% Guardar figura
output_file = 'dispersion_relation.png';
saveas(fig, output_file);
fprintf('Figura guardada: %s\n\n', output_file);


%% Resumen
fprintf('Sin difusión: Re(λ)<0 → estado estable (u0=%.1f, v0=%.1f).\n', u0, v0);
fprintf('Con difusión: lambda_max=%.4f>0 → inestabilidad de Turing.\n', lambda_max);
fprintf('k* = %.4f → L* = %.2f celdas (escala del patrón).\n', k_max, longitud_onda);
fprintf('lambda* debería coincidir con la pendiente lambda del ajuste en Ej.2.\n');

% La inestabilidad de Turing es paradójica: la difusión (normalmente estabilizadora)
% desestabiliza el estado homogéneo para modos con k ≈ k*.
% Mecanismo: u se autocataliza localmente, v suprime a u a distancia (Dv >> Du)
% → manchas de activador con longitud de onda dominante L* = 2π/k*.
