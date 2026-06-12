function [U,epsi,epsi1,epsi3,sigma,sigma1,sigma3] = solid3D(nodes,elements,loads,bcs)

    %% SETUP

    % Number of nodes
    n_nodes = size(nodes, 1);

    % Number of elements
    n_elements = size(elements, 1);
    

    %% ASSEMBLY

    % K-Matrix: Global linear stiffness matrix
    % No dense matrix allocation for efficiency purposes later on
    K_i = zeros(n_elements*144, 1); % row indices of each nonzero entry
    K_j = zeros(n_elements*144, 1); % column indices of each nonzero entry
    K_v = zeros(n_elements*144, 1); % the actual values at those positions
    
    for i = 1:n_elements
        % Nodes: access current element nodes
        node_ids = elements(i, 2:5);        % [n1 n2 n3 n4] as a vector
        node_coords = nodes(node_ids, 2:4); % coords(1,:) = [x1 y1 z1], etc.
                                            % 4x3 matrix, rows are nodes
        
        % Material constants for current element
        E = elements(i,6);
        nu = elements(i,7); % Poisson's ratio
        
        % K-Matrix Assembly
        [Be, De, Ve] = element_matrices(node_coords, E, nu);
        ke = Ve * Be' * De * Be;
        
        dofs = reshape([3*node_ids-2; 3*node_ids-1; 3*node_ids], 1, []);
        
        idx = (i-1)*144 + 1 : i*144;
        [di, dj] = meshgrid(dofs, dofs);
        K_i(idx) = di(:);
        K_j(idx) = dj(:);
        K_v(idx) = ke(:);
    end
    

    %% F_ext: LOAD VECTOR
    F = zeros(3*n_nodes, 1);
    
    for i = 1:size(loads, 1)
        node_id = loads(i, 1);
        dofs    = 3*node_id-2 : 3*node_id;   % [x y z] DOFs for this node
        F(dofs) = loads(i, 2:4);
    end


    %% BOUNDARY CONDITIONS
    n_constrained = sum(~isnan(bcs(:,2:4)), 'all');
    constrained = zeros(n_constrained, 1);
    c_idx = 1;

    for i = 1:size(bcs,1)
        node_id = bcs(i,1);
        r    = bcs(i,2:4);
        dofs = [3*node_id-2, 3*node_id-1, 3*node_id];

        for j = 1:3
            if ~isnan(r(j))
                d = dofs(j);

                % F correction - find triplets in column d
                col_mask = (K_j == d);
                F(K_i(col_mask)) = F(K_i(col_mask)) - K_v(col_mask) * r(j);

                % zero row d - find triplets where row is d
                K_v(K_i == d) = 0;

                % zero column d - find triplets where col is d
                K_v(K_j == d) = 0;

                % set prescribed value
                F(d) = r(j);
                
                % for additional diagonal triplet values
                constrained(c_idx) = d;
                c_idx = c_idx + 1;
            end
        end
    end

    % add diagonal entries all at once
    K_i = [K_i; constrained];
    K_j = [K_j; constrained];
    K_v = [K_v; ones(numel(constrained),1)];

    % K-Matrix: Build sparse K
    K = sparse(K_i, K_j, K_v, 3*n_nodes, 3*n_nodes);


    %% SOLVE
    U = K\F;


    %% POSTPROCESSING

    % Stress and strain allocations
    epsi = zeros(n_elements,7);
    epsi(:,1) = 1:n_elements;
    sigma = zeros(n_elements,7);
    sigma(:,1) = 1:n_elements;

    epsi1  = zeros(n_elements, 1);
    epsi3  = zeros(n_elements, 1);
    sigma1 = zeros(n_elements, 1);
    sigma3 = zeros(n_elements, 1);


    for i = 1:n_elements
        % Nodes: access current element nodes
        node_ids = elements(i, 2:5);        
        node_coords = nodes(node_ids, 2:4);
        dofs = reshape([3*node_ids-2; 3*node_ids-1; 3*node_ids], 1, []);
        Ue = U(dofs);
        
        % Material constants for current element
        E = elements(i,6);
        nu = elements(i,7); 
        
        % B- and D-Matrix recomputation 
        [Be, De,  ~] = element_matrices(node_coords, E, nu);
        
        % Strain and Stress calulation
        eps = Be * Ue;
        sig = De * eps;

        epsi(i,2:7) = eps';
        sigma(i,2:7) = sig';

        % Principal strains
        ex  = epsi(i,2); ey  = epsi(i,3); ez  = epsi(i,4);
        gyz = epsi(i,5); gxz = epsi(i,6); gxy = epsi(i,7);

        % Strain Tensor
        strain_tensor = [ex      gxy/2   gxz/2;
                    gxy/2   ey      gyz/2;
                    gxz/2   gyz/2   ez   ];

        % Eigenvalue extraction
        ev = sort(eig(strain_tensor));
        epsi1(i) = ev(3);   % max principal strain
        epsi3(i) = ev(1);   % min principal strain

        % Principal stresses
        sx  = sigma(i,2); sy  = sigma(i,3); sz  = sigma(i,4);
        tyz = sigma(i,5); txz = sigma(i,6); txy = sigma(i,7);
    
        % Stress tensor
        stress_tensor = [sx  txy txz;
                    txy sy  tyz;
                    txz tyz sz];
    
        % Eigenvalue extraction
        sv = sort(eig(stress_tensor));
        sigma1(i) = sv(3);  % max principal stress
        sigma3(i) = sv(1);  % min principal stress
    end
end

%% SUBFUNCTIONS
function [Be, De, Ve] = element_matrices(node_coords, E, nu)
    % Jacobian: relates physical and reference coordinate derivatives (Kochmann §15.2)
    % Je maps physical->reference, Ji=inv(Je) maps reference->physical
    dx = node_coords(2:4,:) - node_coords(1,:);  % 3x3, differences from node 1
    Je = dx';
    Ve  = det(Je) / 6;

    dN_ref = [ 1  0  0 -1;
               0  1  0 -1;
               0  0  1 -1];   % 3x4, col a -> node a (reference for Na calc)

    dN_phys = Je \ dN_ref;    % 3x4, inv(Je) * dN_ref
    
    % B-Matrix: holds derivatives of Na
    Be = zeros(6,12);
    for a = 1:4
        col = (a-1)*3 + 1;
        dNx = dN_phys(1,a);
        dNy = dN_phys(2,a);
        dNz = dN_phys(3,a);
        Be(:, col:col+2) = [dNx  0    0 ; % εx  = dNx * ux
                           0    dNy  0  ; % εy  = dNy * uy
                           0    0    dNz; % εz  = dNz * uz
                           0    dNz  dNy; % γyz = dNy*uz + dNz*uy 
                           dNz  0    dNx; % γxz = dNx*uz + dNz*ux
                           dNy  dNx  0 ]; % γxy = dNx*uy + dNy*ux
    end

    % D-Matrix: 3D isotropic elasticity stiffness matrix
    De = E/((1+nu)*(1-2*nu)) * ...
        [1-nu  nu    nu    0           0           0; 
         nu    1-nu  nu    0           0           0;
         nu    nu    1-nu  0           0           0;
         0     0     0     (1-2*nu)/2  0           0;
         0     0     0     0           (1-2*nu)/2  0;
         0     0     0     0           0           (1-2*nu)/2];
end