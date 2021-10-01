close all;
clear all;
clc;

% Adaptive mesh refinement using MARIE electromagnetic solver

% The adaptive mesh refinment is performed similarly to HFSS, where the
% mesh is refined at a particular frequency until convergence; the criteria
% for convergence is the maximum change in the S-parameters from iteration
% to iteration

% by Eugene Milshteyn
% Jan 2019

%Coarse mesh to start iteration process
filename = '2port_rect_loop_w_refined_shield'; 

%Filename for refined meshes (will be appended with 2, 3, etc.)
filename2 = '2port_rect_loop_rect_shield_abs_delta_J_length_constraint';

% -------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------
iteration = 1;

max_num_iterations = 15; %number of iterations

freq = 125e6; % 3T; frequency to do adaptive meshing at

Sparam = zeros(2,2,max_num_iterations);

%how many faces of mesh to refine each iteration
percent_refine(1,1) = 25;
percent_refine(1,2:max_num_iterations) = 25;

delta_S_store = zeros(1,max_num_iterations);
delta_S = 1;

delta_S_tol = 0.005; %tolerance of delta S between two iterations

% Refined coil used for snapping coil mesh every iteration
[COIL2] = Import_COIL('2port_rect_loop_rect_shield_high_res.smm');
vertices_STL = COIL2.node;
faces_STL = COIL2.elem;

% -------------------------------------------------------------------------
% Load the corresponding body model
% -------------------------------------------------------------------------

% we use the same way of loading as standard GUI

RHBM_name = 'Sphere_R10cm_5mm_res';

RHBMfile = strcat('/autofs/cluster/wald_deepbrain/eugene/MARIE-master/data/bodies/',[RHBM_name],'.vmm');
    

% -------------------------------------------------------------------------
% Set resolution
% -------------------------------------------------------------------------

resolution = 5e-3; % 5mm

% we just set some dimensions to limit the size of the model
[RHBM] = Import_RHBM(RHBMfile,resolution);

while delta_S > delta_S_tol
    
    disp(['iteration = ', num2str(iteration)])
    
    if iteration < 2
        [COIL] = Import_COIL([filename,'.smm']); %Import the coil for MR_solver
    else
        [COIL] = Import_COIL([filename2,'_',num2str(iteration),'.smm']);
    end
    

    % -------------------------------------------------------------------------
    %                   initialize stuff
    % -------------------------------------------------------------------------

    % generate EM constants
    EMconstants;

    % initialize figure counter
    figidx = 1;

    % find scatterer 
    [L,M,N,~] = size(RHBM.r);
    idxS = find(abs(RHBM.epsilon_r(:)-1 + RHBM.sigma_e(:)) > 1e-12); % non-air positions
    nS = length(idxS);
    nD = L*M*N; % number of voxels in the complete domain
    idxSS = [idxS; nD+idxS; 2*nD+idxS]; % the vector of non-air positions in 3D grid

    % coordinates of domain
    xd = RHBM.r(:,:,:,1);
    yd = RHBM.r(:,:,:,2);
    zd = RHBM.r(:,:,:,3);
    Dcoord = [xd(:), yd(:), zd(:)];

    % scatterer
    xs = xd(idxS);
    ys = yd(idxS);
    zs = zd(idxS);
    Scoord = [xs(:), ys(:), zs(:)];

    figure(figidx);
    plot3(xs, ys, zs ,'bs');
    hold off;
    axis equal;

    % -------------------------------------------------------------------------
    %                 define EM vars and constants
    % -------------------------------------------------------------------------

    % properties at the given frequency
    e_r = RHBM.epsilon_r - 1j*RHBM.sigma_e/(eo*omega);

    % to simplify
    r = RHBM.r; 

    % get the voxel side
    dx = r(2,1,1,1) - r(1,1,1,1);
    % Gram matrix (value)
    Gram = dx^3; % volume of the voxel

    % -------------------------------------------------------------------------
    %                  Generate circulants
    % -------------------------------------------------------------------------

    % compute the circulants
    [fN] = getOPERATORS(r,freq,'N',[],'DEMCEM');
    [fK] = getOPERATORS(r,freq,'K',[],'DEMCEM');
    
    % Flag 6 is for calculating B fields, Flag 9 is for including air voxels
    flags(1:9) = 1; 
    flags(9) = 0;
    
    % call the MARIE solver with the coil
    tol = 1e-7;
    [ZPw,Jcw,Jbw,Sbw,Ebw,Bbw,Gsarw,Pabsw] = MR_Solver(RHBM,COIL,freq,tol,flags);

%     fprintf(1, '\n          Elapsed time  = %g [sec]', toc(tini));
%     fprintf(1, '\n');
    Sparam_temp = z2s(ZPw);
    Sparam(:,:,iteration) = Sparam_temp;
    
    % -------------------------------------------------------------------------
    % Done
    % -------------------------------------------------------------------------

    % the domain r is a L x M x N x 3
    %   with L x M x N the voxels, and the 3 coordinates per voxel
    %
    % idxS gives the position of the body voxels in the domain
    %
    % ZPw is the port Z-parameter of the coils (Nchannels x Nchannels x Num_freqs)
    % Jcw is the coil current functions (Nsegments x Nchannels x Num_freqs)
    % Jbw is the volumetric current distribution in the body (L x M x N x 3 x Nchannels x Num_freqs)
    %   it has  3 components per voxel
    %   and as many vectors as number of channels
    %   and number of analysis (or frequencies, in this case Num_freqs = 1
    % Sbw is the SAR (L x M x N x Nchannels x Num_freqs)
    % Ebw is the E field distribution in the domain (L x M x N x 3 x Nchannels x Num_freqs)
    % Bbw is the B field distribution in the domain (L x M x N x 3 x Nchannels x Num_freqs)
    % Gsarw is the global SAR (Nchannels x Num_freqs)
    % Pabsw is the absorbed Power (Nchannels x Num_freqs)
    %
    % For more details see the MARIE documentation
    %

    if iteration < 2
        SOLfile = strcat([filename],'_',[RHBM_name], '_solution.mat');
    else
        SOLfile = strcat([filename2],'_',num2str(iteration),'_',[RHBM_name], '_solution.mat');
    end
    % -------------------------------------------------------------------------
    %                  Calculate delta_S
    % -------------------------------------------------------------------------
    disp('Calculating delta_S')
    
    if iteration == 1
        delta_S = 1;
    else
        delta_S = max(max(abs(Sparam(:,:,iteration)-Sparam(:,:,iteration-1))));
    end
      
    delta_S_store(1,iteration) = delta_S;
      
    if delta_S < delta_S_tol
        continue;
    end
    % -------------------------------------------------------------------------
    %                  Select faces to refine using Surface Currents
    % -------------------------------------------------------------------------
    
    Jcoil = Jcw;
    
    Ct = COIL.Ct;
    
    [j_cart_ct, Element_areas] = calculate_surface_currents_V3(COIL,Jcoil);

    % Calculating total current at each face
    j_combined = j_cart_ct(:,:,1); %single port
    j_total = abs(sqrt(j_combined(:,1).^2+j_combined(:,2).^2+j_combined(:,3).^2)).*Element_areas';
    
    % Created interpolation function for previous coil
    if iteration > 1
        F = scatteredInterpolant(Ct_previous(:,1),Ct_previous(:,2),Ct_previous(:,3),j_total_previous);
        temp = Ct';
        interpolated_current = F(temp(:,1),temp(:,2),temp(:,3));
        delta_J = abs(interpolated_current-j_total);
        face_surf_curr = [(1:size(Ct,2))' delta_J];
    else
        face_surf_curr = [(1:size(Ct,2))' j_total];
        delta_J = 0;
    end
    
    j_total_previous = j_total;
    Ct_previous = Ct';
    
    % Sort faces and currents from highest current to lowest current
    faces_resorted = sortrows(face_surf_curr,2,'descend');
    
    % Calculate edge lengths within each face and find smallest edge length
    for ii = 1:size(faces_resorted,1)
        p1 = COIL.node(:,COIL.elem(1,ii));
        p2 = COIL.node(:,COIL.elem(2,ii));
        p3 = COIL.node(:,COIL.elem(3,ii));
        dist = [norm(p1-p2) norm(p1-p3) norm(p2-p3)];
        faces_resorted(ii,3) = min(dist)*1000; %convert to mm from m
    end
    
    ind = find(faces_resorted(:,3)>=1.25);
    
    faces_too_small = find(faces_resorted(:,3)<1.25);
    
    faces_resorted_2 = faces_resorted(ind,:);
    
    % Select percentage of faces to refine, starting with ones with highest
    % gradient
    
    num_faces_to_refine = ceil(size(faces_resorted_2,1)*percent_refine(1,iteration)/100);
    
    % Select the faces that see the highest E gradient; choose top X%
    
    faces_to_refine = faces_resorted_2(1:num_faces_to_refine,1);
    
    faces_not_to_refine = [faces_resorted((num_faces_to_refine+1):size(faces_resorted_2,1),1) faces_too_small];  
    
    
    
    save(SOLfile, 'COIL', 'RHBM', 'freq', 'ZPw', 'Jcw', 'Jbw', 'Sbw', 'Ebw', 'Bbw', 'Gsarw', 'Pabsw','faces_to_refine', 'faces_not_to_refine', 'delta_S_store','delta_J','Sparam','Sparam_temp', '-v7.3');
      
      
    % -------------------------------------------------------------------------
    %                   Refine the mesh
    % -------------------------------------------------------------------------
    disp('Refining mesh')
    
    if iteration == 1
        refine_MARIE_mesh_V11(filename,filename2,iteration+1,faces_to_refine,faces_not_to_refine,vertices_STL,faces_STL);
    else
        refine_MARIE_mesh_V11([filename2,'_',num2str(iteration)],filename2,iteration+1,faces_to_refine,faces_not_to_refine,vertices_STL,faces_STL);
    end
    
    iteration = iteration + 1;

    % Do not exceed max_num_iterations
    if iteration > max_num_iterations
        continue;
    end
    
    
end

