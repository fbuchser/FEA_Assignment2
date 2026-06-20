%% SETTINGS
OutputPath = 'results_ansys\compare\';
if ~exist('results_ansys\', 'dir'), mkdir('results_ansys'); end
if ~exist(OutputPath, 'dir'), mkdir('results_ansys\compare'); end

% Set Latex design
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(groot, 'defaultColorbarTickLabelInterpreter', 'latex')
set(groot, 'defaultLegendInterpreter', 'latex')

% Custom colour map
n = 256;
r = [linspace(1,1,n/2), linspace(1,0,n/2)]';
g = [linspace(0,1,n/2), linspace(1,0,n/2)]';
b = [linspace(0,1,n/2), linspace(1,1,n/2)]';
coolwarm = [r, g, b];
coolwarm = flipud(coolwarm);


%% DATA PREP - Custom solver
disp('Loading Data');
load('data\FemurInput.mat')           % nodes, elements
load('results\sideways_fall_results.mat')  % U, sigma1, sigma3, epsi1, epsi3, bcs, loads

n_nodes = size(nodes, 1);
n_elem  = size(elements, 1);

node_ids_solver = nodes(:,1);
elem_ids_solver = elements(:,1);

conn = elements(:,2:5);
coords = nodes(:,2:4);

U_mat = reshape(U, 3, n_nodes)';

%% PARSE ANSYS FILES
disp('Parsing ANSYS files');
[ansys_disp_ids,   ansys_disp_vals]   = parse_nodal_listing(  'results_ansys\numerical\displacements.lis', 3);
[ansys_stress_ids, ansys_stress_vals] = parse_element_listing('results_ansys\numerical\prin_stress.lis',   n_elem);
[ansys_strain_ids, ansys_strain_vals] = parse_element_listing('results_ansys\numerical\prin_strain.lis',   n_elem);
%% ALIGN TO SOLVER ORDERING
U_ansys = zeros(n_nodes, 3);
[~, ia, ib] = intersect(node_ids_solver, ansys_disp_ids);
U_ansys(ia,:) = ansys_disp_vals(ib,:);

s1_ansys = NaN(n_elem,1); [~,ia,ib] = intersect(elem_ids_solver, ansys_stress_ids); s1_ansys(ia) = ansys_stress_vals(ib,1);
s3_ansys = NaN(n_elem,1); [~,ia,ib] = intersect(elem_ids_solver, ansys_stress_ids); s3_ansys(ia) = ansys_stress_vals(ib,3);
e1_ansys = NaN(n_elem,1); [~,ia,ib] = intersect(elem_ids_solver, ansys_strain_ids); e1_ansys(ia) = ansys_strain_vals(ib,1);
e3_ansys = NaN(n_elem,1); [~,ia,ib] = intersect(elem_ids_solver, ansys_strain_ids); e3_ansys(ia) = ansys_strain_vals(ib,3);

%% RELATIVE ERROR FIELDS
disp('Computing relative error fields');
u_mag_own   = sqrt(sum(U_mat.^2,   2));
u_mag_ansys = sqrt(sum(U_ansys.^2, 2));

bc_nodes = unique(bcs(:,1));
bc_node_idx = ismember(node_ids_solver, bc_nodes);

tmp = (u_mag_own - u_mag_ansys) ./ u_mag_ansys * 100;
tmp(bc_node_idx) = NaN;
du_mag_elem = mean(tmp(conn), 2, 'omitnan');
nan_mask = isnan(du_mag_elem);
if any(nan_mask)
    neighbor_vals = mean(du_mag_elem(conn(nan_mask,:)), 2, 'omitnan');
    du_mag_elem(nan_mask) = neighbor_vals;
end

tmp = U_mat(:,3) - U_ansys(:,3);
duz_elem = mean(tmp(conn), 2);

ds1 = sigma1 - s1_ansys;
ds3 = sigma3 - s3_ansys;
de1 = epsi1  - e1_ansys;
de3 = epsi3  - e3_ansys;


%% PLOTS

 areplot_diff(du_mag_elem, conn, coords, coolwarm, OutputPath, 'Relative error $[\%]$', 'Displacement magnitude')
plot_diff(duz_elem,    conn, coords, coolwarm, OutputPath, 'Absolute error $[mm]$', 'Z displacement')
plot_diff(ds1,         conn, coords, coolwarm, OutputPath, 'Absolute error $[N/mm^2]$', 'Max principal stress')
plot_diff(ds3,         conn, coords, coolwarm, OutputPath, 'Absolute error $[N/mm^2]$', 'Min principal stress')
plot_diff(de1,         conn, coords, coolwarm, OutputPath, 'Absolute error', 'Max principal strain')
plot_diff(de3,         conn, coords, coolwarm, OutputPath, 'Absolute error', 'Min principal strain')


%% HELPERS

function plot_diff(data, conn, coords, coolwarm, OutputPath, mode_str, title_str)

disp(['Plotting ' title_str]);
t_plot = tic;


fig = figure('Name', title_str, 'Color', 'w', 'Visible', 'off');

t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'loose');
t.OuterPosition = [0 0 1 0.95];

nexttile
tetramesh(conn, coords, data); colormap(coolwarm); 
axis equal; view(-150, 30); % FRONT
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k')
t_title = title('Anterior'); t_title.FontSize = 9;
t_title.Units = 'normalized'; t_title.Position(2) = 0.9;
xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]')

nexttile
tetramesh(conn, coords, data); colormap(coolwarm);
axis equal; view(75, 30); % BACK
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k')
t_title = title('Lateral-Posterior'); t_title.FontSize = 9;
t_title.Units = 'normalized'; t_title.Position(2) = 0.95;
xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]')

node_acc = accumarray(conn(:), repmat(data, 4, 1), [size(coords,1) 1], @sum, 0);
node_cnt = accumarray(conn(:), ones(4*size(conn,1), 1), [size(coords,1) 1], @sum, 0);
node_vals = node_acc ./ max(node_cnt, 1);
clim_val = max(abs(node_vals(isfinite(node_vals))));
set(findall(fig, 'type', 'axes'), 'CLim', [-clim_val, clim_val]);
cb = colorbar; cb.Color = 'k';
sgtitle([mode_str ': ' title_str ' - Sideways Fall'], 'Color', 'k', 'interpreter','latex', 'FontWeight', 'normal')

fprintf('Plot time: %.2f s\n', toc(t_plot))
t_plot = tic;

disp(['Saving ' title_str]);

file_str = strrep(title_str, ' ', '_');
exportgraphics(fig, fullfile(OutputPath, [file_str '.png']), 'BackgroundColor', 'w', 'Resolution', 300)
close(fig)

fprintf('Save time: %.2f s\n', toc(t_plot))
end


function [node_ids, vals] = parse_nodal_listing(filepath, n_cols)
fid = fopen(filepath, 'r');
node_ids = zeros(5000, 1);
vals     = zeros(5000, n_cols);
count    = 0;
while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line), break; end
    nums = regexp(strtrim(line), '[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?', 'match');
    if numel(nums) < 1 + n_cols, continue; end
    id = str2double(nums{1});
    if isnan(id) || id ~= round(id) || id <= 0, continue; end
    row = str2double(nums(2:1+n_cols));
    if any(isnan(row)), continue; end
    count = count + 1;
    node_ids(count)  = id;
    vals(count,:)    = row;
end
fclose(fid);
node_ids = node_ids(1:count);
vals     = vals(1:count,:);
end


function [elem_ids, vals] = parse_element_listing(filepath, n_elem)
fid = fopen(filepath, 'r');
elem_ids     = zeros(n_elem, 1);
vals         = zeros(n_elem, 3);
current_rows = zeros(10, 3);
current_id   = [];
n_rows       = 0;
count        = 0;

    function flush()
        if ~isempty(current_id) && n_rows > 0
            count = count + 1;
            elem_ids(count) = current_id;
            vals(count,:)   = mean(current_rows(1:n_rows,:), 1);
        end
    end

while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line), break; end
    tok = regexp(line, 'ELEMENT=\s*(\d+)', 'tokens');
    if ~isempty(tok)
        flush();
        current_id = str2double(tok{1}{1});
        n_rows     = 0;
        continue
    end
    nums = regexp(strtrim(line), '[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?', 'match');
    if numel(nums) < 4, continue; end
    id = str2double(nums{1});
    if isnan(id) || id ~= round(id) || id <= 0, continue; end
    row = str2double(nums(2:4));
    if any(isnan(row)), continue; end
    n_rows = n_rows + 1;
    current_rows(n_rows,:) = row;
end
flush();
fclose(fid);
elem_ids = elem_ids(1:count);
vals     = vals(1:count,:);
end

