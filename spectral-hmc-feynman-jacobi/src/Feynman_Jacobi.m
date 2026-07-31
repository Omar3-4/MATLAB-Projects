function [Classical_Path, Quantum_Fluctuations, L_class, L_quantum_min] = Feynman_Jacobi(xa,xb,T,m,Vx,dVx,D,N,L,NO_Test,Plot_Flag,Pot_Name)
% FEYNMAN_JACOBI Analyzes classical paths and quantum fluctuations via HMC
%
% Inputs:
%   xa, xb      - Boundary conditions (start and end points)
%   T           - Total time
%   m           - Mass
%   Vx, dVx     - Potential function handle and its gradient
%   D           - Number of time slices
%   N           - Number of HMC samples
%   L           - Number of leapfrog steps for HMC
%   NO_Test     - Number of test trajectories for Jacobi metric
%   Plot_Flag   - Boolean flag to generate visualizations
%   Pot_Name    - String name of the potential for plot titles

    % PROTECTION: Ensure required inputs exist to prevent nargin crashes
    if nargin < 10
        error('Feynman_Jacobi requires at least 10 inputs.');
    end
    if nargin < 11 || isempty(Plot_Flag), Plot_Flag = false; end
    if nargin < 12 || isempty(Pot_Name), Pot_Name = 'Custom Potential'; end
    
    % Ensure utils/Jacobi_metric is on MATLAB path
    if exist('Jacobi_metric', 'file') ~= 2
        utils_dir = fullfile(fileparts(mfilename('fullpath')), 'utils');
        if exist(fullfile(utils_dir, 'Jacobi_metric.m'), 'file')
            addpath(utils_dir);
        end
    end

    try
        % PROTECTION: Catch setup errors or missing HMC function dependency
        dt = T/(D-1);
        x_base_full = linspace(xa, xb, D);
        x_base_mid = x_base_full(2:end-1);
        Vx_shifted = @(x_f) Vx(x_f + x_base_mid(:));
        dVx_shifted = @(x_f) dVx(x_f + x_base_mid(:));
        [Samples_f,~,~] = HMC(Vx_shifted,dVx_shifted,D-2,N,L,m,dt);
    catch ME
        % FALLBACK: Return NaN arrays if setup or HMC crashes
        Classical_Path = NaN(1, D); Quantum_Fluctuations = NaN(1, D);
        L_class = NaN; L_quantum_min = NaN; return;
    end
    
    try
        % PROTECTION: Catch matrix dimension mismatches during path construction
        zeros_vec = zeros(size(Samples_f,1), 1);
        Samples_fluct_full = [zeros_vec, Samples_f, zeros_vec];
        x = Samples_fluct_full + repmat(x_base_full, size(Samples_f,1), 1);
        
        Classical_Path = mean(x, 1);
        Quantum_Fluctuations = std(x, 0, 1);
        
        x_mid = (Classical_Path(1:end-1) + Classical_Path(2:end)) / 2;
        K_Ec = 0.5*m*(diff(Classical_Path)/dt).^2;
        P_Ec = Vx(x_mid);
        E_array = K_Ec + P_Ec;
        E_true = mean(E_array);
        
        L_class = Jacobi_metric(Classical_Path,E_true,m,Vx);
    catch ME
        % FALLBACK: Abort computation if statistical processing fails
        Classical_Path = NaN(1, D); Quantum_Fluctuations = NaN(1, D);
        L_class = NaN; L_quantum_min = NaN; return;
    end
    
    try
        % PROTECTION: Catch out-of-bounds indices in random path sampling
        x_q = x(randperm(size(x,1),NO_Test),:);
        L_q = zeros(1,NO_Test);
        for i = 1:NO_Test
            L_q(i) = Jacobi_metric(x_q(i,:),E_true,m,Vx);
        end
        L_quantum_min = min(L_q);
    catch
        % FALLBACK: Skip quantum metric bounds if path extraction fails
        L_q = NaN(1, NO_Test);
        L_quantum_min = NaN;
    end

    fprintf('Length of Classical Geodesic : %.4f\n', L_class);
    fprintf('Minimum Quantum Path Length  : %.4f\n', L_quantum_min);
    
    if (nargin >= 11) && Plot_Flag == true
        if nargin < 12
            Pot_Name = 'Custom Potential';
        end
        
        try
            % PROTECTION: Isolate graphics processing to prevent UI runtime crashes
            x_axis = 1:D;
            figure('Name', sprintf('Feynman-Jacobi HMC Analysis [%s]', Pot_Name), 'NumberTitle', 'off', 'Units', 'normalized','Position', [0.2, 0.2, 0.6, 0.6]);
            
            subplot(2, 3, 1);
            hold on; grid on;
            fill([x_axis, fliplr(x_axis)], ...
                 [Classical_Path + Quantum_Fluctuations, fliplr(Classical_Path - Quantum_Fluctuations)], ...
                 [0.8 0.8 0.8], 'EdgeColor', 'none');
            plot(x_axis, Classical_Path, 'r', 'LineWidth', 2);
            title('Quantum vs Classical Dynamics');
            xlabel('Time Slices'); ylabel('Position x(\tau)');
            
            subplot(2, 3, 2);
            random_idx = randi(size(x,1));  
            plot(x_axis, x(random_idx, :), 'b', 'LineWidth', 1);
            title('Single Quantum Trajectory (Instantons)');
            xlabel('Time Slices'); ylabel('Position');
            grid on;
            
            subplot(2, 3, 3);
            mid_point_idx = round(D/2);
            plot(1:size(x,1), x(:, mid_point_idx), 'k', 'LineWidth', 0.5);
            title(sprintf('HMC Trace Plot (at t = %d)', mid_point_idx));
            xlabel('HMC Iteration'); ylabel('Position');
            grid on;
            
            subplot(2, 3, 4);
            plot(1:D-1, E_array, 'b', 'LineWidth', 1.5);
            yline(E_true, 'r--', 'LineWidth', 2);
            title(sprintf('Conservation of Energy (Mean = %.4f)', E_true));
            xlabel('Time Slices'); ylabel('Total Energy E');
            axis tight; grid on;
            
            subplot(2, 3, 5);
            X_mat = repmat(x_axis(:), size(x,1), 1);
            Y_mat = x(:);
            hist3([X_mat, Y_mat], [100, 100], 'CDataMode','auto','FaceColor','interp', 'EdgeAlpha', 0);
            view(2); 
            colormap jet;
            hold on;
            plot3(x_axis, Classical_Path, max(histcounts(Y_mat))*ones(1,D), 'r', 'LineWidth', 2);
            title('Quantum Path Density Heatmap');
            xlabel('Time Slices'); ylabel('Position');
            
            subplot(2, 3, 6);
            hold on; grid on;
            histogram(L_q, 'FaceColor', [0.6 0.6 0.6]);
            xline(L_class, 'r', 'LineWidth', 3);
            title('Jacobi Metric: Absolute Minimum Proof');
            xlabel('Jacobi Arc Length (L)'); ylabel('Frequency');
            
            sgtitle(sprintf('Comprehensive HMC Path Integral Analysis\nPotential: %s', Pot_Name), 'FontSize', 16, 'FontWeight', 'bold');
            
            figure('Name', sprintf('Potential Landscape & Stats [%s]', Pot_Name), 'NumberTitle', 'off', 'Units', 'normalized','Position', [0.2, 0.2, 0.6, 0.6]);
            x_min_plot = min(min(x)) - 1;
            x_max_plot = max(max(x)) + 1;
            x_grid = linspace(x_min_plot, x_max_plot, 1000);
            V_grid = Vx(x_grid);
            
            subplot(1, 3, 1);
            plot(x_grid, V_grid, 'k', 'LineWidth', 2);
            hold on; grid on;
            plot(xa, Vx(xa), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
            plot(xb, Vx(xb), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
            xline(min(Classical_Path), 'b--');
            xline(max(Classical_Path), 'b--');
            title('The Potential Landscape V(x)');
            xlabel('Position x'); ylabel('Potential Energy V(x)');
            legend('V(x)', 'Start (A)', 'End (B)', 'Path Min/Max', 'Location', 'best');
            
            subplot(1, 3, 2);
            [T_mesh, X_mesh] = meshgrid(linspace(0, T, 50), linspace(x_min_plot, x_max_plot, 100));
            V_mesh = Vx(X_mesh);
            surf(T_mesh, X_mesh, V_mesh, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
            colormap parula;
            hold on;
            time_actual = linspace(0, T, D);
            V_path = Vx(Classical_Path);
            plot3(time_actual, Classical_Path, V_path, 'r', 'LineWidth', 3);
            view(-30, 30);
            title('3D Classical Trajectory on V(x)');
            xlabel('Time'); ylabel('Position x'); zlabel('V(x)');
            grid on;
            
            subplot(1, 3, 3);
            mid_pt_trace = x(:, round(D/2)); 
            N_lags = min(100, length(mid_pt_trace)-1);
            acf = zeros(N_lags+1, 1);
            trace_mean = mean(mid_pt_trace);
            trace_var = var(mid_pt_trace);
            for lag = 0:N_lags
                cov_val = mean((mid_pt_trace(1:end-lag) - trace_mean) .* (mid_pt_trace(lag+1:end) - trace_mean));
                acf(lag+1) = cov_val / trace_var;
            end
            stem(0:N_lags, acf, 'Marker', 'none', 'LineWidth', 1.5);
            yline(0, 'k-');
            yline(0.05, 'r--'); yline(-0.05, 'r--');
            title('HMC Autocorrelation (Middle Point)');
            xlabel('Lag (Iterations)'); ylabel('ACF');
            grid on;
            
            
        catch
            % FALLBACK: Ignore plotting errors to ensure calculation returns cleanly
            warning('Plotting phase encountered an error and was bypassed.');
        end
    end
end