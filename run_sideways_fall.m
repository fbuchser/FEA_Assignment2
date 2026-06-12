%% LOAD INPUT

load('data\FemurInput.mat')


%% BOUNDARY CONDITIONS

% Shaft support
MinShaftZ = min(nodes(:,4));
I = nodes(:,4) < MinShaftZ+5;
bc_nodes1 = nodes(I,1);
nn = length(bc_nodes1);
bc1 = [bc_nodes1 NaN(nn,1) NaN(nn,1) zeros(nn,1)];
    
% Head support
MinHeadX = min(nodes(:,2));
I = nodes(:,2) < MinHeadX+5;
bc_nodes2 = nodes(I,1);
nn = length(bc_nodes2);
bc2 = [bc_nodes2 zeros(nn,1) NaN(nn,1) NaN(nn,1)];
    
% Greater trochanter support
MaxGTX = max(nodes(:,2));
I = nodes(:,2) > MaxGTX-5;
bc_nodes3 = nodes(I,1);
nn = length(bc_nodes3);
bc3 = [bc_nodes3 -1*ones(nn,1) zeros(nn,1) NaN(nn,1)];

bcs = [bc1; bc2; bc3]; % -1mm: displacement of 1 mm at the greater trochanter 

%% LOADS
loads = [1 0 0 0]; % no external force — loading is via prescribed displacement at GT


%% RUN SOLVER
[U, epsi, epsi1, epsi3, sigma, sigma1, sigma3] = solid3D(nodes, elements, loads, bcs);


%% SAVE OUTPUT
OutputPath = 'results\';
if ~exist(OutputPath, 'dir'), mkdir('results'); end

save(fullfile(OutputPath, 'sideways_fall_results.mat'), 'U', 'epsi', 'epsi1', 'epsi3', 'sigma', 'sigma1', 'sigma3')


%% PLOT SOLVER RESULTS
plot_femur_results(nodes, elements, U, epsi1, epsi3, OutputPath, 'Sideways Fall')