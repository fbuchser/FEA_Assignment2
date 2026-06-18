%% LOAD INPUT
disp('Loading Femur Input');

load('data\FemurInput.mat')


%% SURFACE NODES
% Surface node extraction

% Extract all 4 faces from each tet
faces = [elements(:,[2,3,4]); elements(:,[2,3,5]); elements(:,[2,4,5]); elements(:,[3,4,5])];

% Sort each face's nodes
faces_sorted = sort(faces, 2);

% Find surface faces (appear exactly once)
[~, ~, ic] = unique(faces_sorted, 'rows');
counts = accumarray(ic, 1);
surface_faces = faces_sorted(counts(ic) == 1, :);

% Surface nodes
surface_nodes = unique(surface_faces(:));
surface_node_coords = nodes(ismember(nodes(:,1), surface_nodes), :);


%% BOUNDARY CONDITIONS
disp('Setting Boundary Conditions');

% Shaft support
MinShaftZ = min(nodes(:,4));
I = surface_node_coords(:,4) < MinShaftZ+0.5;
bc_nodes1 = surface_node_coords(I,1);
nn = length(bc_nodes1);
bcs = [bc_nodes1 zeros(nn,1) zeros(nn,1) zeros(nn,1)];


%% LOADS
disp('Setting Loads');

% - Force
team_weight =  65 + 59;
F_total = team_weight * 9.81;

% 7° from shaft axis in frontal plane, 28° in transverse plane
% (Bergmann et al. 2001, Fig. 8, standing on one leg)
Fz =  F_total * cosd(7);             % proximal (along shaft)
Fx = -F_total * sind(7) * cosd(28);  % medial = negative X
Fy = -F_total * sind(7) * sind(28);  % anterior direction = negative Y

% - Node Selection

% Medio-cranial reference points
[~, I_minX] = min(nodes(:,2));
[~, I_maxZ] = max(nodes(:,4));
X_of_maxZ = nodes(I_maxZ, 2);
Z_of_minX = nodes(I_minX, 4);

% Filter surface nodes by medio-cranial criteria
I = surface_node_coords(:,2) < X_of_maxZ & surface_node_coords(:,4) > Z_of_minX;

% Limit Y-range
Y_coords = surface_node_coords(I, 3);
Y_mid = (max(Y_coords) + min(Y_coords)) / 2;
Y_range = (max(Y_coords) - min(Y_coords)) * 0.66;
I = I & surface_node_coords(:,3) > Y_mid - Y_range*0.5 & surface_node_coords(:,3) < Y_mid + Y_range*0.5;

head_nodes = surface_node_coords(I, 1);
nn_head = length(head_nodes);

loads = [head_nodes, (Fx/nn_head)*ones(nn_head,1), (Fy/nn_head)*ones(nn_head,1), (Fz/nn_head)*ones(nn_head,1)];


%% RUN SOLVER
disp('Running Solver');

[U, epsi, epsi1, epsi3, sigma, sigma1, sigma3] = solid3D(nodes, elements, loads, bcs);


%% SAVE OUTPUT
disp('Saving Solver Results');

OutputPath = 'results\';
if ~exist(OutputPath, 'dir'), mkdir('results'); end

save(fullfile(OutputPath, 'single_leg_stance_results.mat'), 'bcs', 'loads', 'U', 'epsi', 'epsi1', 'epsi3', 'sigma', 'sigma1', 'sigma3')


%% PLOT SOLVER RESULTS
disp('-- Plotting Solver Results --');

plot_femur_results(nodes, elements, U, sigma1, sigma3, epsi1, epsi3, bcs, loads, OutputPath, 'Single Leg Stance')
