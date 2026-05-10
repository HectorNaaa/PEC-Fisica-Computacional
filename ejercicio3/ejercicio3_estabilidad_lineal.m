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

% =========================================================================
%% SECCIÓN 0 — PARÁMETROS DEL MODELO
% =========================================================================

a  = 0.1;     % Parámetro de alimentación del Brusselator
b  = 0.9;     % Parámetro de drenaje del Brusselator
Du = 1.0;     % Coeficiente de difusión del activador u (lento)
Dv = 10.0;    % Coeficiente de difusión del inhibidor v (rápido, Dv >> Du)

% ---- Estado estacionario homogéneo ----
% Condición: f(u0,v0) = 0  y  g(u0,v0) = 0
%   a - u0 + u0²·v0 = 0
%   b - u0²·v0      = 0  →  u0²·v0 = b
% Sumando: a - u0 + b = 0  →  u0 = a + b
%          v0 = b / u0²
u0 = a + b;
v0 = b / u0^2;

fprintf('=================================================================\n');
fprintf('  EJERCICIO 3 — Análisis de estabilidad lineal (Brusselator)\n');
fprintf('=================================================================\n');
fprintf('Parámetros:\n');
fprintf('  a  = %.2f,  b  = %.2f\n', a, b);
fprintf('  Du = %.2f,  Dv = %.2f\n', Du, Dv);
fprintf('\nEstado estacionario homogéneo:\n');
fprintf('  u0 = a + b       = %.6f\n', u0);
fprintf('  v0 = b / (a+b)²  = %.6f\n', v0);
fprintf('=================================================================\n\n');


% =========================================================================
%% SECCIÓN 1 — JACOBIANO DEL SISTEMA SIN DIFUSIÓN
% =========================================================================
%
% La estabilidad sin difusión se determina linealizando el sistema ODE
% alrededor de (u0, v0). Las derivadas parciales de f y g son:
%
%   ∂f/∂u = -1 + 2·u·v    ∂f/∂v =  u²
%   ∂g/∂u = -2·u·v        ∂g/∂v = -u²
%
% El Jacobiano J = [∂f/∂u  ∂f/∂v ; ∂g/∂u  ∂g/∂v] evaluado en (u0,v0):

J11 = -1 + 2*u0*v0;   % ∂f/∂u = b - a     (para estos valores = 0.8)
J12 =  u0^2;           % ∂f/∂v = (a+b)²   (= 1.0)
J21 = -2*u0*v0;        % ∂g/∂u = -(b-a)   (= -0.8... espera, = -2*1*0.9 = -1.8... hmm)
J22 = -u0^2;           % ∂g/∂v = -(a+b)²  (= -1.0)

% Nota: con a=0.1, b=0.9: u0=1, v0=0.9
%   J11 = -1 + 2·1·0.9 = 0.8
%   J12 =  1² = 1
%   J21 = -2·1·0.9 = -1.8
%   J22 = -1² = -1

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


% =========================================================================
%% SECCIÓN 2 — AUTOVALORES SIN DIFUSIÓN
% =========================================================================
%
% Los autovalores λ de J determinan cómo evolucionan perturbaciones
% homogéneas (k=0):
%
%   perturbación ~ exp(λ·t)
%
%   Re(λ) < 0 para ambos λ → la perturbación decae → ESTADO ESTABLE
%   Re(λ) > 0 para algún λ → la perturbación crece → ESTADO INESTABLE
%
% Si los autovalores son complejos conjugados α ± βi con α < 0,
% el punto fijo es un FOCO ESTABLE: las perturbaciones oscilan mientras
% decaen (la solución "gira" en espiral hacia el punto fijo).

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

% ─────────────────────────────────────────────────────────────────────────
% Interpretación clave:
% Si el sistema SIN difusión es estable (Re(λ) < 0) pero CON difusión
% aparece inestabilidad, entonces la difusión es responsable del patrón.
% Esto es precisamente la INESTABILIDAD DE TURING o "diffusion-driven
% instability": un paradoja en que el mecanismo estabilizador (difusión)
% desestabiliza el estado homogéneo.
% ─────────────────────────────────────────────────────────────────────────
fprintf('\n');


% =========================================================================
%% SECCIÓN 3 — RELACIÓN DE DISPERSIÓN CON DIFUSIÓN
% =========================================================================
%
% Al añadir difusión, una perturbación espacial de la forma:
%
%   δu(x,t) ~ exp(λ·t) · cos(k·x)
%
% evoluciona según la matriz modificada por difusión:
%
%   M(k) = J - k²·D    donde  D = diag(Du, Dv)
%
%         = [ J11 - Du·k²    J12        ]
%           [ J21            J22 - Dv·k² ]
%
% El número de onda k (en rad/unidad de longitud) está relacionado con
% la longitud de onda espacial λ_onda por: λ_onda = 2π/k
%
% Para cada k calculamos el autovalor dominante (mayor parte real).
% Si ese autovalor tiene Re > 0, el modo k crece exponencialmente.

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


% =========================================================================
%% SECCIÓN 4 — MÁXIMO DE λ(k): MODO DOMINANTE Y ESCALA ESPACIAL
% =========================================================================
%
% El máximo de λ(k) sobre todos los k es el modo que crece más rápido.
% Este modo domina el patrón emergente y determina su escala espacial.
%
%   k* = argmax λ(k)              → número de onda dominante
%   L* = 2π / k*                  → longitud de onda dominante
%
% Si k* es grande → estructuras pequeñas (manchas finas)
% Si k* es pequeño → estructuras grandes (manchas anchas)

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


% =========================================================================
%% SECCIÓN 5 — FIGURA: RELACIÓN DE DISPERSIÓN
% =========================================================================
%
% La relación de dispersión λ(k) muestra:
%   • k < k1 o k > k2 → λ < 0 → modos estables (decaen)
%   • k1 < k < k2     → λ > 0 → modos inestables (crecen) ← BANDA DE TURING
%   • k = k*          → λ máxima → modo dominante

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


% =========================================================================
%% SECCIÓN 6 — RESUMEN FINAL
% =========================================================================
fprintf('=================================================================\n');
fprintf('  RESUMEN — Resultados esperados\n');
fprintf('=================================================================\n');
fprintf('\n');
fprintf('1. SIN DIFUSIÓN:\n');
fprintf('   El estado estacionario (u0=%.1f, v0=%.1f) es estable.\n', u0, v0);
fprintf('   Re(λ) < 0 para ambos autovalores de J.\n');
fprintf('   Las perturbaciones homogéneas decaen → sin patrón posible.\n');
fprintf('\n');
fprintf('2. CON DIFUSIÓN:\n');
fprintf('   Aparece una banda de k con λ(k) > 0 (inestabilidad de Turing).\n');
fprintf('   Esto ocurre porque Dv >> Du: el inhibidor v difunde mucho más\n');
fprintf('   rápido que el activador u, rompiendo el balance espacial.\n');
fprintf('\n');
fprintf('3. MODO DOMINANTE:\n');
fprintf('   k* = %.4f rad/píxel  →  L* = 2π/k* ≈ %.2f píxeles.\n', k_max, longitud_onda);
fprintf('   El patrón de manchas en la simulación del Ej.1 debería tener\n');
fprintf('   una separación media entre manchas ≈ L* = %.2f celdas.\n', longitud_onda);
fprintf('\n');
fprintf('4. CONEXIÓN CON EL EJERCICIO 1:\n');
fprintf('   La tasa λ* predicha aquí debe coincidir con la pendiente λ\n');
fprintf('   medida empíricamente en el ajuste de A(t) del Ejercicio 2.\n');
fprintf('\n');
fprintf('=================================================================\n');

% ─────────────────────────────────────────────────────────────────────────
% NOTA PARA LA DEFENSA ORAL:
%
% La inestabilidad de Turing es contraintuitiva: la difusión, que
% normalmente ESTABILIZA, aquí DESESTABILIZA el estado homogéneo para
% modos espaciales concretos. El mecanismo es:
%   1. El activador u se autocataliza localmente (término u²v).
%   2. El inhibidor v difunde rápido y suprime u a distancia.
%   3. Localmente, u "gana" antes de ser suprimido → manchas de activador.
%   4. A distancia, v "gana" → fondo con poco activador.
% El resultado es un patrón espacial estacionario con longitud de onda L*.
% ─────────────────────────────────────────────────────────────────────────
