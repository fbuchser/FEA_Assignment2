%% LOAD INPUT
disp('Loading Femur Input');

load('data\FemurInput.mat')


%% BOUNDARY CONDITIONS
disp('Setting Boundary Conditions');

% Shaft support
MinShaftZ = min(nodes(:,4));
I = nodes(:,4) < MinShaftZ+5;
bc_nodes1 = nodes(I,1);
nn = length(bc_nodes1);
bcs = [bc_nodes1 zeros(nn,1) zeros(nn,1) zeros(nn,1)];


%% LOADS
disp('Setting Loads');

team_weight =  65 + 59;
F_total = team_weight * 9.81;

% 7° from shaft axis in frontal plane, 28° in transverse plane
% (Bergmann et al. 2001, Fig. 8, standing on one leg)
Fz =  F_total * cosd(7);             % proximal (along shaft)
Fx = -F_total * sind(7) * cosd(28);  % medial = negative X
Fy = -F_total * sind(7) * sind(28);  % anterior direction = negative Y

MaxHeadZ = max(nodes(:,4));
I_head = nodes(:,4) < MaxHeadZ+5;
head_nodes = nodes(I_head,1);
nn_head = length(head_nodes);

loads = [head_nodes, (Fx/nn_head)*ones(nn_head,1), (Fy/nn_head)*ones(nn_head,1), (Fz/nn_head)*ones(nn_head,1)];

%% RUN SOLVER
disp('Running Solver');

[U, epsi, epsi1, epsi3, sigma, sigma1, sigma3] = solid3D(nodes, elements, loads, bcs);


%% SAVE OUTPUT
disp('Saving Solver Results');

OutputPath = 'results\';
if ~exist(OutputPath, 'dir'), mkdir('results'); end

save(fullfile(OutputPath, 'single_leg_stance_results.mat'), 'bcs', 'U', 'epsi', 'epsi1', 'epsi3', 'sigma', 'sigma1', 'sigma3')


%% PLOT SOLVER RESULTS
disp('-- Plotting Solver Results --');

plot_femur_results(nodes, elements, U, epsi1, epsi3, bcs, OutputPath, 'Single Leg Stance')
