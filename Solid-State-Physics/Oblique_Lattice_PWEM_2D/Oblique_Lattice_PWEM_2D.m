% =========================================================================
% Method: 2D Plane Wave Expansion Method (PWEM)
% Model: Parallelogram / Square (at a=0) Lattice Geometry 
% Outputs:
%   1. 2D Electronic Band Structure along Gamma-X-M-Y-Gamma
%   2. Real-Space Parallelogram Unit Cell Geometry & Primitive Vectors
%   3. 1st Brillouin Zone & Reciprocal Cell with High-Symmetry Points
%   4. Constant Energy Surface (CES) (Band 1)
% =========================================================================

tic
clear; clc; close all;

%% 1. Geometric & Physical Parameters 
l = 0.645;             % length in nm
a = 0.200;             % shift in nm
h = 0.554;             % height in nm
b = l - a;             % b in nm
c = a;
V0 = 2.0;              % Potential barrier in eV
omega = l * h;
d = sqrt(a^2 + h^2);

% example: Effective mass for ReS2: m* = 0.45 * m_e
% hm = hbar^2 / (2 * m*) in (eV * nm^2)
hm = 0.0381 / 0.45;    % ~ 0.0847 eV*nm^2

%% 2. Plane Wave Basis Setup
N = 10;                     % Maximum plane wave index
m1 = -N:N; 
m2 = -N:N;
M_len = length(m1);              
Mk = (2*N+1)^2;                % Number of k-points for plotting
k = linspace(-pi/d, pi/d, Mk); % Wave vector array across the 1st Brillouin Zone
[M1, M2] = meshgrid(m1, m2);
m1_all = M1(:);
m2_all = M2(:);
Gx = (2*pi/omega)*m1_all* h;
Gy = (2*pi/omega)*(m2_all*l - m1_all*a);
N_pw = length(m1_all);

% Reciprocal basis vectors
b1 = [ (2*pi)/l ; -(2*pi*a)/(l*h) ];
b2 = [ 0        ; (2*pi)/h        ];

% High Symmetry points
Gamma = [0; 0];
X     = 0.5 * b1;
M     = 0.5 * (b1 + b2);
Y     = 0.5 * b2;

sym_points = [Gamma, X, M, Y, Gamma]; 
num_segments = size(sym_points, 2) - 1;
N_interp = 50;

% Linear Interpolation on Gamma-X-M-Y-Gamma path
k_path = [];
for seg = 1:num_segments
    k_start = sym_points(:, seg);
    k_end   = sym_points(:, seg+1);
    t = linspace(0, 1, N_interp);
    seg_k = (1 - t) .* k_start + t .* k_end;
    if seg == 1
        k_path = [k_path, seg_k];
    else
        k_path = [k_path, seg_k(:, 2:end)]; 
    end
end

N_total = size(k_path, 2);
d_path = zeros(1, N_total);
for j = 2:N_total
    delta_k = k_path(:, j) - k_path(:, j-1);
    d_path(j) = d_path(j-1) + norm(delta_k);
end

d_ticks = zeros(1, num_segments + 1);
d_ticks(1) = 0;

current_idx = 1;
for seg = 1:num_segments
    current_idx = current_idx + (N_interp - 1);
    d_ticks(seg + 1) = d_path(current_idx);
end

% 3. Fourier Components of Potential V
V = zeros(N_pw);
for ii = 1:N_pw
    for j = 1:N_pw
        dm1 = m1_all(ii)-m1_all(j); % Difference between plane wave indices
        dm2 = m2_all(ii)-m2_all(j);
        if dm1==0 && dm2==0
            V(ii,j) = (2*V0/omega) * (l+d);
        elseif dm1==0 && dm2~=0
            V(ii,j) = (2*V0/omega) * l;
        elseif dm1~=0 && dm2==0
            V(ii,j) = (2*V0/omega) * d;
        elseif dm1~=0 && dm2~=0
            V(ii,j) = 0;
        end
    end
end
V = (V + V') / 2; % Enforce exact Hermiticity

% 4. Solve Band Structure along k-path
num_bands = 2;
E_all = zeros(num_bands, N_total); 

for j = 1:N_total
    KE = diag(hm * ((k_path(1, j) + Gx).^2 + (k_path(2, j) + Gy).^2));
    H = KE + V;
    H = (H + H') / 2;
    eig_vals = sort(real(eig(H)));
    E_all(:, j) = eig_vals(1:num_bands);
end

%% Plotting

% Figure 1: Electronic Band Structure with Shaded Gap
figure('Color', 'w', 'Name', 'Electronic_Band_Structure');
hold on;

% Shading band gap between Band 1 and Band 2
top_band1 = max(E_all(1, :));
bot_band2 = min(E_all(2, :));
if bot_band2 > top_band1
    fill([0, d_ticks(end), d_ticks(end), 0], ...
         [top_band1, top_band1, bot_band2, bot_band2], ...
         [0.92, 0.92, 0.92], 'EdgeColor', 'none', 'DisplayName', 'Band Gap');
end

colors = lines(num_bands);
for b_idx = 1:num_bands
    plot(d_path, E_all(b_idx, :), 'Color', colors(b_idx, :), 'LineWidth', 2, ...
         'DisplayName', sprintf('Band %d', b_idx));
end

for i = 1:length(d_ticks)
    xline(d_ticks(i), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
end

xlim([0, d_ticks(end)]);
set(gca, 'XTick', d_ticks, 'XTickLabel', {'\Gamma', 'X', 'M', 'Y', '\Gamma'}, 'FontSize', 11);
ylabel('Energy (eV)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('k-path', 'FontSize', 12, 'FontWeight', 'bold');
title('Electronic Band Structure along High-Symmetry Path', 'FontSize', 13);
legend('Location', 'northeast');
grid on; box on;


% Figure 2: Real Space Unit Cell Geometry (Dynamic & Exact Vector Polygon)
poly_x = [0, l, l + a, a, 0];
poly_y = [0, 0, h, h, 0];

figure('Color', 'w', 'Name', 'Real_Space_Unit_Cell');
hold on;
fill(poly_x, poly_y, [0.85, 0.92, 1], 'EdgeColor', 'b', 'LineWidth', 2, 'DisplayName', 'Unit Cell');
plot(poly_x, poly_y, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5, 'HandleVisibility', 'off');

% Real-space lattice vectors (a1 along x, a2 to (a,h))
quiver(0, 0, l, 0, 0, 'k', 'LineWidth', 1.8, 'MaxHeadSize', 0.2, 'DisplayName', 'a_1');
quiver(0, 0, a, h, 0, 'k', 'LineWidth', 1.8, 'MaxHeadSize', 0.2, 'DisplayName', 'a_2');

span_x = max(poly_x) - min(poly_x);
span_y = max(poly_y) - min(poly_y);
pad_x = 0.15 * span_x;
pad_y = 0.15 * span_y;

text(l/2, -0.06*span_y, '\bf a_1', 'FontSize', 11, 'HorizontalAlignment', 'center');
text(a/2 - 0.04*span_x, h/2 + 0.04*span_y, '\bf a_2', 'FontSize', 11, 'HorizontalAlignment', 'center');

xlim([min(poly_x) - pad_x, max(poly_x) + pad_x]);
ylim([min(poly_y) - pad_y, max(poly_y) + pad_y]);
axis equal; grid on; box on;
if a==0
    title('Real Space Unit Cell (Square)', 'FontSize', 13);
else 
    title('Real Space Unit Cell (Parallelogram)', 'FontSize', 13)
end
xlabel('x (nm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y (nm)', 'FontSize', 12, 'FontWeight', 'bold');


% Figure 3: Reciprocal Unit Cell & 1st Brillouin Zone with High-Symmetry Path
cell_x = [0, b1(1), b1(1)+b2(1), b2(1), 0];
cell_y = [0, b1(2), b1(2)+b2(2), b2(2), 0];

bz_x = [M(1), M(1)-b1(1), -M(1), -M(1)+b1(1), M(1)];
bz_y = [M(2), M(2)-b1(2), -M(2), -M(2)+b1(2), M(2)];

figure('Color', 'w', 'Name', 'Brillouin_Zone_and_Reciprocal_Cell');
hold on;

fill(cell_x, cell_y, [1 0.9 0.9], 'EdgeColor', 'r', 'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Reciprocal Cell');
fill(bz_x, bz_y, [0.85 0.92 1], 'EdgeColor', 'b', 'LineWidth', 2, 'DisplayName', '1st BZ');

% Trace k-path inside BZ
plot(k_path(1, :), k_path(2, :), 'Color', [0.8 0 0], 'LineWidth', 2, 'DisplayName', 'k-path');

quiver(0, 0, b1(1), b1(2), 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.2, 'HandleVisibility', 'off');
quiver(0, 0, b2(1), b2(2), 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.2, 'HandleVisibility', 'off');

plot(0, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'HandleVisibility', 'off');
plot(X(1), X(2), 'ms', 'MarkerFaceColor', 'm', 'MarkerSize', 6, 'HandleVisibility', 'off');
plot(M(1), M(2), 'gd', 'MarkerFaceColor', 'g', 'MarkerSize', 6, 'HandleVisibility', 'off');
plot(Y(1), Y(2), 'c^', 'MarkerFaceColor', 'c', 'MarkerSize', 6, 'HandleVisibility', 'off');

all_kx = [cell_x, bz_x];
all_ky = [cell_y, bz_y];
span_kx = max(all_kx) - min(all_kx);
span_ky = max(all_ky) - min(all_ky);
dx = 0.04 * span_kx;
dy = 0.04 * span_ky;

text(b1(1)/2 + dx, b1(2)/2 - dy, '\bf b_1', 'FontSize', 12);
text(b2(1)/2 - 1.5*dx, b2(2)/2 + dy, '\bf b_2', 'FontSize', 12);
text(0 + dx, 0 - dy, '\bf \Gamma', 'FontSize', 13);
text(X(1) + dx, X(2), '\bf X', 'FontSize', 11, 'Color', 'm');
text(M(1) + dx, M(2) + dy, '\bf M', 'FontSize', 11, 'Color', 'g');
text(Y(1) - 1.5*dx, Y(2) + dy, '\bf Y', 'FontSize', 11, 'Color', 'c');

pad_kx = 0.15 * span_kx;
pad_ky = 0.15 * span_ky;
xlim([min(all_kx) - pad_kx, max(all_kx) + pad_kx]);
ylim([min(all_ky) - pad_ky, max(all_ky) + pad_ky]);
axis equal;

xlabel('k_x (nm^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('k_y (nm^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
title('1st Brillouin Zone & Reciprocal Cell with k-Path', 'FontSize', 13);
legend('Location', 'northwest');
grid on; box on;


% Figure 4: Constant Energy Surface (CES) / Iso-energy contours (Proportional Scaling)
Nk_2d = 60;
span_bz_x = max(bz_x) - min(bz_x);
span_bz_y = max(bz_y) - min(bz_y);
pad_ces_x = 0.15 * span_bz_x;
pad_ces_y = 0.15 * span_bz_y;

kx_vec = linspace(min(bz_x) - pad_ces_x, max(bz_x) + pad_ces_x, Nk_2d);
ky_vec = linspace(min(bz_y) - pad_ces_y, max(bz_y) + pad_ces_y, Nk_2d);
[KX, KY] = meshgrid(kx_vec, ky_vec);

E_grid = zeros(Nk_2d, Nk_2d);
target_band = 1; % 1st Band

for iy = 1:Nk_2d
    for ix = 1:Nk_2d
        kx_i = KX(iy, ix);
        ky_i = KY(iy, ix);
        KE_i = diag(hm * ((kx_i + Gx).^2 + (ky_i + Gy).^2));
        H_i = KE_i + V;
        eig_vals_i = sort(real(eig(H_i)));
        E_grid(iy, ix) = eig_vals_i(target_band);
    end
end

figure('Color', 'w');
hold on;
contourf(KX, KY, E_grid, 25, 'LineColor', 'none');
colorbar;
colormap(jet);

plot(bz_x, bz_y, 'w--', 'LineWidth', 2);
plot(0, 0, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 6);

xlim([min(kx_vec), max(kx_vec)]);
ylim([min(ky_vec), max(ky_vec)]);
axis equal;
xlabel('k_x (nm^{-1})'); ylabel('k_y (nm^{-1})');
title(['Constant Energy Surface (Band ' num2str(target_band) ')']);
grid on; box on;


%% Save Figures 
[script_folder, script_name] = fileparts(mfilename('fullpath'));
if isempty(script_name)
    script_name = 'Parallelogram_PWEM';
    script_folder = pwd;
end
output_dir = fullfile(script_folder, script_name);

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fig_handles = findobj('Type', 'figure');
for i_fig = 1:length(fig_handles)
    fig = fig_handles(i_fig);
    if isvalid(fig)
        fig_name = get(fig, 'Name');
        if isempty(fig_name)
            fig_name = sprintf('figure_%d', get(fig, 'Number'));
        else
            fig_name = regexprep(fig_name, '[^\w\-_]', '_');
        end
        file_path = fullfile(output_dir, sprintf('%s.png', fig_name));
        exportgraphics(fig, file_path, 'Resolution', 300);
    end
end

toc
