%%% Simulación principal para el Modelo de Crecimiento de Placa de Aterosclerosis
%%% Adaptado para usabilidad y reproducibilidad en GNU Octave
%%% Basado en: Friedman, Hao, & Hu (2015)
%%%========================================================================

function test3()
    % 1. Parámetros (Uso de estructura para evitar variables globales)
    p.k_1 = 1.4;        % Tasa de ingestión de LDL por macrófagos (/día)
    p.k_2 = 10;         % Tasa de consumo de HDL (/día)
    p.K_1 = 1e-2;       % Constante de Michaelis-Menten para LDL (g/cm^3)
    p.K_2 = 0.5;        % Constante de Michaelis-Menten para HDL (g/cm^3)
    p.mu_1 = 0.003;     % Tasa de muerte de macrófagos (/día)
    p.mu_2 = 0.005;     % Tasa de muerte de células espumosas (/día)
    p.r_1 = 2.4e-5;     % Tasa de degradación de LDL (/día)
    p.r_2 = 5.5e-7;     % Tasa de degradación de HDL (/día)
    p.lambda = 2.57e-3; % Tasa de activación de macrófagos (día)
    p.delta = 2.54e-5;  % Factor de inhibición (g/cm^3)
    p.M_0 = 5e-4;       % Densidad basal de macrófagos (g/cm^3)
    p.L_0 = 200e-5;     % Concentración de LDL en sangre (g/cm^3)
    p.H_0 = 40e-5;      % Concentración de HDL en sangre (g/cm^3)

    % 2. Condiciones Iniciales [L, H, M, F]
    z_ini = [p.L_0, p.H_0, p.M_0, 0]; 
    tspan = [0, 300];           % Simulación de 300 días

    % 3. Resolver las EDOs
    % Uso de función anónima para pasar la estructura de parámetros.
    % Los "Function handles" (@) son la sintaxis preferida en GNU Octave y MATLAB.
    try
        [t, z] = ode15s(@(t, z) odefun(t, z, p), tspan, z_ini);
    catch
        % Respaldo para entornos donde ode15s no está disponible o falla
        warning('ode15s falló, intentando con ode45...');
        [t, z] = ode45(@(t, z) odefun(t, z, p), tspan, z_ini);
    end

    % 4. Procesamiento de Resultados
    % R representa la carga de placa normalizada (M+F relativo al M0 basal)
    R = (z(:,3) + z(:,4)) ./ p.M_0;

    % 5. Graficación
    figure(1);
    % Verificar existencia de sgtitle (MATLAB R2018b+ u Octave reciente)
    if exist('sgtitle', 'file') || exist('sgtitle', 'builtin')
        sgtitle('Evolución Temporal de los Componentes de la Placa');
    end

    titles = {'Concentración de LDL', 'Concentración de HDL', ...
              'Densidad de Macrófagos', 'Densidad de Células Espumosas'};
    ylabels = {'L (g/cm^3)', 'H (g/cm^3)', 'M (g/cm^3)', 'F (g/cm^3)'};

    for i = 1:4
        subplot(2,2,i);
        plot(t, z(:,i), 'LineWidth', 1.5);
        title(titles{i}); xlabel('Tiempo (días)'); ylabel(ylabels{i});
        grid on;
    end

    figure(2);
    plot(t, R, 'g', 'LineWidth', 2);
    xlabel('Tiempo (días)'); ylabel('Carga de Placa Normalizada R(t)');
    title(sprintf('Dinámica de Estabilidad: L_0 = %.2e, H_0 = %.2e', p.L_0, p.H_0));
    grid on;
end

function dz = odefun(t, z, p)
    % Desempaquetar variables de estado
    L = z(1); H = z(2); M = z(3); F = z(4);
    dz = zeros(4,1);

    % Ecuaciones Diferenciales (Friedman et al., 2015)
    dz(1) = p.L_0 - p.k_1*M*L/(p.K_1+L) - p.r_1*L; % LDL
    dz(2) = p.H_0 - p.k_2*H*F/(p.K_2+F) - p.r_2*H; % HDL
    dz(3) = - p.k_1*M*L/(p.K_1+L) + p.k_2*H*F/(p.K_2+F) + ...
             p.lambda*M*L/(H-p.delta) - p.mu_1*M; % Macrófagos
    dz(4) = p.k_1*M*L/(p.K_1+L) - p.k_2*H*F/(p.K_2+F) - p.mu_2*F; % Células Espumosas
end