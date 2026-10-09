%%% =========================================================================
%%% Graficador de Zonas de Riesgo de Aterosclerosis (Sección 12.1)
%%% Adaptado para máxima compatibilidad con GNU Octave y MATLAB
%%% Basado en: Friedman, Hao, & Hu (2015) / Chou & Friedman (2016)
%%% Evalúa la razón de placa R(300) = (M(300) + F(300)) / M_0 sobre una malla de
%%% concentraciones sanguíneas iniciales de LDL (L_0) y HDL (H_0).
%%% =========================================================================

function risk()
close all;

% Detección de entorno GNU Octave y selección de toolkit estable
isOctave = (exist('OCTAVE_VERSION', 'builtin') ~= 0);
if isOctave
    avail_tk = available_graphics_toolkits();
    if any(strcmp('qt', avail_tk))
        graphics_toolkit('qt');
    elseif any(strcmp('gnuplot', avail_tk))
        graphics_toolkit('gnuplot');
    end
end

%% 1. Definición de Parámetros (Estructura de parámetros portable)
p.k1     = 1.4;                 % Tasa de ingestión de LDL por macrófagos (/día)
p.k2     = 10;                  % Tasa de remoción de LDL oxidado por HDL (/día)
p.K1     = 1e-2;                % Constante de semisaturación (g/cm^3)
p.K2     = 0.5;                 % Constante de semisaturación (g/cm^3)
p.mu1    = 0.003;               % Tasa de muerte de macrófagos (/día)
p.mu2    = 0.005;               % Tasa de muerte de células espumosas (/día)
p.r1     = 2.4e-5;              % Tasa de degradación de LDL (/día)
p.r2     = 5.5e-7;              % Tasa de degradación de HDL (/día)
p.lambda = 2.57e-3;             % Coeficiente de infiltración (/día)
p.delta  = 2.54e-5;             % Umbral de inhibición (g/cm^3)
p.M0     = 5e-4;                % Densidad inicial de macrófagos (g/cm^3)

T_final  = 300;                 % Período de tiempo (días)

%% 2. Definir la Malla para L_0 y H_0 (en unidades de 10^-5 g/cm^3)
L0_vec = linspace(50, 250, 40);  % Concentraciones de LDL
H0_vec = linspace(20, 100, 40);  % Concentraciones de HDL

[L0_grid, H0_grid] = meshgrid(L0_vec, H0_vec);
R_final = zeros(size(L0_grid));

%% 3. Integración Numérica sobre la Malla
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% Selección automática de solver compatible con Octave / MATLAB
if exist('ode15s', 'file') || exist('ode15s', 'builtin')
    ode_solver = @ode15s;
else
    ode_solver = @ode45;
end

fprintf('Calculando la razón de carga de placa R(300) sobre la malla (L0, H0)...\n');
total_rows = size(L0_grid, 1);
for i = 1:total_rows
    % Notificar progreso y ceder tiempo a la interfaz gráfica para evitar bloqueos
    if mod(i, 5) == 0 || i == 1 || i == total_rows
        fprintf('  Progreso: fila %d/%d (%.1f%%)...\n', i, total_rows, (i/total_rows)*100);
        drawnow;
    end
    for j = 1:size(L0_grid, 2)
        % Convertir valores de la malla a concentraciones reales (g/cm^3)
        L0_val = L0_grid(i,j) * 1e-5;
        H0_val = H0_grid(i,j) * 1e-5;
        
        % Condiciones iniciales: [L(0), H(0), M(0), F(0)]
        z_ini = [L0_val, H0_val, p.M0, 0];
        
        % Resolver el sistema de EDOs
        try
            [~, z] = ode_solver(@(t, z) fun_atherosclerosis(t, z, L0_val, H0_val, p), ...
                [0, T_final], z_ini, options);
        catch
            [~, z] = ode45(@(t, z) fun_atherosclerosis(t, z, L0_val, H0_val, p), ...
                [0, T_final], z_ini, options);
        end
        
        % Calcular R(300) = w(300) / M_0
        w_final = z(end, 3) + z(end, 4);
        R_final(i, j) = w_final / p.M0;
    end
end
fprintf('Cálculo finalizado exitosamente.\n');

%% 4. Graficación de las Zonas de Riesgo y Libres de Riesgo
h_fig = figure(1);
set(h_fig, 'Name', 'Zona de Riesgo de Aterosclerosis', 'Color', 'w');

% Mapa de contornos de los valores de R(300)
[~, h_c] = contourf(L0_grid, H0_grid, R_final, 25, 'LineStyle', 'none');
hold on;
cb = colorbar;
try
    cb.Label.String = 'Razón Relativa de Carga de Placa R(300)';
catch
    ylabel(cb, 'Razón Relativa de Carga de Placa R(300)');
end
colormap(jet);

% Línea de contorno que separa las zonas de Riesgo (R > 1) y Libre de Riesgo (R < 1)
[C, h_b] = contour(L0_grid, H0_grid, R_final, [1 1], 'w--', 'LineWidth', 2.5);
try
    clabel(C, h_b, 'FontSize', 10, 'Color', 'white', 'FontWeight', 'bold');
catch
    try
        clabel(C, h_b);
    catch
        % Fallback silencioso si el toolkit gráfico no soporta clabel
    end
end

% Graficar Puntos de Referencia AHA mencionados en la Sección 12
plot(190, 40, 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
text(193, 38, 'Alto Riesgo (190, 40)', 'Color', 'white', 'FontWeight', 'bold');

plot(110, 50, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
text(113, 52, 'Libre de Riesgo (110, 50)', 'Color', 'white', 'FontWeight', 'bold');

xlabel('L_0 (Nivel de LDL en 10^{-5} g/cm^3)', 'FontSize', 12);
ylabel('H_0 (Nivel de HDL en 10^{-5} g/cm^3)', 'FontSize', 12);
title('Mapa de Zonas de Riesgo de Aterosclerosis en T = 300 Días', 'FontSize', 13);
grid on;

legend(h_b, {'Umbral de Riesgo (R = 1)'}, 'Location', 'northwest');

%% 5. Guardado de Figuras Compatible con GNU Octave y MATLAB
folderName = 'fig';
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

png_path = fullfile(folderName, sprintf('fig1_%s.png', timestamp));
eps_path = fullfile(folderName, sprintf('fig1_%s.eps', timestamp));
fig_path = fullfile(folderName, sprintf('fig1_%s.fig', timestamp));

if ~isOctave && (exist('exportgraphics', 'file') || exist('exportgraphics', 'builtin'))
    exportgraphics(gcf, png_path, 'Resolution', 300);
    exportgraphics(gcf, eps_path);
else
    % Exportación de PNG de alta resolución
    try
        print(gcf, png_path, '-dpng', '-r300');
    catch
        try
            print(gcf, png_path, '-dpng');
        catch err
            warning('No se pudo guardar PNG: %s', err.message);
        end
    end

    % Salvaguarda para EPS en GNU Octave:
    % En Octave con toolkit 'fltk', GL2PS sufre de bucles infinitos de búfer con 'contourf'.
    if isOctave
        cur_tk = graphics_toolkit();
        if strcmp(cur_tk, 'fltk')
            warning(['Se omite la exportación EPS en Octave bajo toolkit FLTK para evitar ' ...
                     'bloqueo de GL2PS con contourf. Use toolkit "qt" o "gnuplot" si requiere EPS.']);
            eps_path = 'omitido (incompatible con FLTK/GL2PS)';
        else
            try
                print(gcf, eps_path, '-depsc');
            catch err
                warning('Octave no pudo exportar EPS: %s', err.message);
                eps_path = 'omitido por error';
            end
        end
    else
        print(gcf, eps_path, '-depsc');
    end
end

if exist('savefig', 'file') || exist('savefig', 'builtin')
    try
        savefig(gcf, fig_path);
    catch
        % Respaldo para Octave si savefig difiere
    end
end
fprintf('Figuras guardadas: PNG [%s], EPS [%s]\n', png_path, eps_path);
end

%% 6. Función Local: Sistema de EDOs de Aterosclerosis
function dz = fun_atherosclerosis(~, z, L0, H0, p)
L = z(1);
H = z(2);
M = z(3);
F = z(4);

dz = zeros(4,1);

% dL/dt
dz(1) = L0 - p.k1 * M * (L / (p.K1 + L)) - p.r1 * L;

% dH/dt
dz(2) = H0 - p.k2 * H * (F / (p.K2 + F)) - p.r2 * H;

% dM/dt
dz(3) = - p.k1 * M * (L / (p.K1 + L)) + p.k2 * H * (F / (p.K2 + F)) ...
    + p.lambda * M * L / (H - p.delta) - p.mu1 * M;

% dF/dt
dz(4) = p.k1 * M * (L / (p.K1 + L)) - p.k2 * H * (F / (p.K2 + F)) - p.mu2 * F;
end