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
bc1 = [bc_nodes1 NaN(nn,1) NaN(nn,1) zeros(nn,1)];
    
% Head support
MinHeadX = min(nodes(:,2));
I = surface_node_coords(:,2) < MinHeadX+5;
bc_nodes2 = surface_node_coords(I,1);
nn = length(bc_nodes2);
bc2 = [bc_nodes2 zeros(nn,1) NaN(nn,1) NaN(nn,1)];
    
% Greater trochanter support
MaxGTX = max(nodes(:,2));
I = surface_node_coords(:,2) > MaxGTX-5;
bc_nodes3 = surface_node_coords(I,1);
nn = length(bc_nodes3);
bc3 = [bc_nodes3 -1*ones(nn,1) zeros(nn,1) NaN(nn,1)];

bcs = [bc1; bc2; bc3]; % -1mm: displacement of 1 mm at the greater trochanter 


%% LOADS
disp('Setting Loads');

loads = [1 0 0 0]; % no external force — loading is via prescribed displacement at GT


%% RUN SOLVER
disp('Running Solver');

[U, epsi, epsi1, epsi3, sigma, sigma1, sigma3] = solid3D(nodes, elements, loads, bcs);


%% SAVE OUTPUT
disp('Saving Solver Results');

OutputPath = 'results\';
if ~exist(OutputPath, 'dir'), mkdir('results'); end

save(fullfile(OutputPath, 'sideways_fall_results.mat'), 'bcs', 'loads', 'U', 'epsi', 'epsi1', 'epsi3', 'sigma', 'sigma1', 'sigma3')


%% PLOT SOLVER RESULTS
disp('-- Plotting Solver Results --');

plot_femur_results(nodes, elements, U, sigma1, sigma3, epsi1, epsi3, bcs, loads, OutputPath, 'Sideways Fall')
