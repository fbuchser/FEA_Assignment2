function plot_femur_results(nodes, elements, U, sigma1, sigma3, bcs, loads, OutputPath, title_str)

set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(groot, 'defaultColorbarTickLabelInterpreter', 'latex')
set(groot, 'defaultLegendInterpreter', 'latex')

if ~exist(OutputPath, 'dir'), mkdir(OutputPath); end
file_str = strrep(title_str, ' ', '_');

% Plot general data prep
coords = nodes(:,2:4);
conn   = elements(:,2:5);

% Plot a data prep
bc_marked = coords(bcs(:,1), :);
non_null_loads = loads(any(loads(:,2:4) ~= 0, 2), :);
loads_marked = coords(non_null_loads(:,1), :);

disp(size(bc_marked))
disp(size(bcs))
disp(size(loads_marked))
disp(size(loads))

% Plot c data prep
n_nodes = size(nodes, 1);
U_mat = reshape(U, 3, n_nodes)';
def_mag = sqrt(sum(U_mat.^2, 2));
def_mag_elem = mean(def_mag(conn), 2);

%% (a) Input mesh
disp('Plotting Input Mesh');

t_plot = tic;
figA = figure('Name', [title_str ' - Mesh'], 'Color', 'w', 'Visible', 'off');

tetramesh(conn, coords, 'FaceColor', [0.4 0.8 0.7], 'FaceAlpha', 0.4);

hold on

scatter3(bc_marked(:,1), bc_marked(:,2), bc_marked(:,3), 'o', 'MarkerEdgeColor', 'b')
scatter3(loads_marked(:,1), loads_marked(:,2), loads_marked(:,3), 'o', 'MarkerEdgeColor', 'r')

hold off
axis equal; view(-150, 30);
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k')

title(['Input mesh with marked constrained nodes - ' title_str], 'Color', 'k')
xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]')

fprintf('Plot (a) time: %.2f s\n', toc(t_plot))
t_plot = tic;

disp('Saving Input Mesh');
exportgraphics(figA, fullfile(OutputPath, [file_str '_in_mesh.png']), 'BackgroundColor', 'w', 'Resolution', 300)
close(figA)

fprintf('Save (a) time: %.2f s\n', toc(t_plot))

if strcmp(title_str, 'Sideways Fall') || ~exist(fullfile(OutputPath, 'Input_Edist.png'), 'file')

    %% (b) Input Young's modulus distribution
    disp('Plotting Input E Distribution');
    
    t_plot = tic;
    figB = figure('Name', [title_str ' - E distribution'], 'Color', 'w', 'Visible', 'off');
    
    t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'loose');
    t.OuterPosition = [0 0 1 0.95];
    
    nexttile
    tetramesh(conn, coords, elements(:,6)); colormap(jet); 
    axis equal; view(-150, 30); % FRONT
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k')
    t_title = title('Anterior'); t_title.FontSize = 9;
    t_title.Units = 'normalized'; t_title.Position(2) = 0.9;
    xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]')
    
    nexttile
    tetramesh(conn, coords, elements(:,6)); colormap(jet);
    axis equal; view(60, 30); % BACK
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k')
    t_title = title('Lateral-Posterior'); t_title.FontSize = 9;
    t_title.Units = 'normalized'; t_title.Position(2) = 0.95;
    xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]')
    
    nexttile
    tetramesh(conn, coords, elements(:,6)); colormap(jet);
    axis equal; view(0, -90); % BELOW
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k')
    t_title = title('Distal'); t_title.FontSize = 9;
    t_title.Units = 'normalized'; t_title.Position(2) = 1.2;
    xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]')
    
    cb = colorbar; cb.Color = 'k';
    sgtitle('Input Young''s modulus $[N/mm^2]$ distribution', 'Color', 'k', 'interpreter','latex', 'FontWeight', 'normal')
    
    fprintf('Plot (b) time: %.2f s\n', toc(t_plot))
    t_plot = tic;
    
    disp('Saving Input E Distribution');
    exportgraphics(figB, fullfile(OutputPath, 'Input_Edist.png'), 'BackgroundColor', 'w', 'Resolution', 300)
    close(figB)
    
    fprintf('Save (b) time: %.2f s\n', toc(t_plot))
end


%% (c) Output 

% deformation magnitude
disp('Plotting Output Deformation Magnitude');
t_plot = tic;

figC1 = figure('Name', [title_str ' - Output: Deformation Magnitude'], 'Color', 'w', 'Visible', 'off');

t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'loose');
t.OuterPosition = [0 0 1 0.95];

nexttile
tetramesh(conn, coords, def_mag_elem); colormap(jet); 
axis equal; view(-150, 30); % FRONT
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k')
t_title = title('Anterior'); t_title.FontSize = 9;
t_title.Units = 'normalized'; t_title.Position(2) = 0.9;
xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]')

nexttile
tetramesh(conn, coords, def_mag_elem); colormap(jet);
axis equal; view(75, 30); % BACK
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k')
t_title = title('Lateral-Posterior'); t_title.FontSize = 9;
t_title.Units = 'normalized'; t_title.Position(2) = 0.95;
xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]')

cb = colorbar; cb.Color = 'k';
sgtitle(['Output: Deformation magnitude $[mm]$ - ' title_str], 'Color', 'k', 'interpreter','latex', 'FontWeight', 'normal')

fprintf('Plot (c1) time: %.2f s\n', toc(t_plot))
t_plot = tic;

disp('Saving Output Deformation Magnitude');
exportgraphics(figC1, fullfile(OutputPath, [file_str '_out_defmag.png']), 'BackgroundColor', 'w', 'Resolution', 300)
close(figC1)

fprintf('Save (c1) time: %.2f s\n', toc(t_plot))


if strcmp(title_str, 'Sideways Fall')
    view2_az = 75;
    view2_el = 30;
    view2_title = 'Lateral-Posterior';
    view2_pos = 0.95;
elseif  strcmp(title_str, 'Single Leg Stance')
    view2_az = -30;
    view2_el = 30;
    view2_title = 'Medial-Posterior';
    view2_pos = 0.9;
end


% sigma1
disp('Plotting Output Sigma1');
t_plot = tic; 

figC2 = figure('Name', [title_str ' - Output: Maximal Principal Stress'], 'Color', 'w', 'Visible', 'off');

t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'loose');
t.OuterPosition = [0 0 1 0.95];

nexttile
tetramesh(conn, coords, sigma1); colormap(jet);
axis equal; view(-150, 30);
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k')
t_title = title('Anterior'); t_title.FontSize = 9;
t_title.Units = 'normalized'; t_title.Position(2) = 0.9;
xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]')

nexttile
tetramesh(conn, coords, sigma1); colormap(jet);
axis equal; view(view2_az, view2_el);
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k')
t_title = title(view2_title); t_title.FontSize = 9;
t_title.Units = 'normalized'; t_title.Position(2) = view2_pos;
xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]')

cb = colorbar; cb.Color = 'k';
sgtitle(['Output: Maximal principal stress $[N/mm^2]$ - ' title_str], 'Color', 'k', 'interpreter','latex', 'FontWeight', 'normal')

fprintf('Plot (c2) time: %.2f s\n', toc(t_plot))
t_plot = tic;

disp('Saving Output Sigma1');
exportgraphics(figC2, fullfile(OutputPath, [file_str '_out_sigma1.png']), 'BackgroundColor', 'w', 'Resolution', 300)
close(figC2)

fprintf('Save (c2) time: %.2f s\n', toc(t_plot))


% sigma3 
disp('Plotting Output Sigma3');
t_plot = tic; 

figC3 = figure('Name', [title_str ' - Output: Minimal Principal Stress'], 'Color', 'w', 'Visible', 'off');

t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'loose');
t.OuterPosition = [0 0 1 0.95];

nexttile
tetramesh(conn, coords, sigma3); colormap(jet);
axis equal; view(-150, 30);
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k')
t_title = title('Anterior'); t_title.FontSize = 9;
t_title.Units = 'normalized'; t_title.Position(2) = 0.9;
xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]')

nexttile
tetramesh(conn, coords, sigma3); colormap(jet);
axis equal; view(view2_az, view2_el);
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k')
t_title = title(view2_title); t_title.FontSize = 9;
t_title.Units = 'normalized'; t_title.Position(2) = view2_pos;
xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]')

cb = colorbar; cb.Color = 'k';
sgtitle(['Output: Minimal principal stress $[N/mm^2]$ - ' title_str], 'Color', 'k', 'interpreter','latex', 'FontWeight', 'normal')

fprintf('Plot (c3) time: %.2f s\n', toc(t_plot))
t_plot = tic;

disp('Saving Output Sigma3');
exportgraphics(figC3, fullfile(OutputPath, [file_str '_out_sigma3.png']), 'BackgroundColor', 'w', 'Resolution', 300)
close(figC3)

fprintf('Save (c3) time: %.2f s\n', toc(t_plot))

end