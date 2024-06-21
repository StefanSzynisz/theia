function [U_return,Epsilon_return,Sigma_return] = forward_problem_2x2_elems_ES(E_arg,nu_arg,P_arg,varargin)
% Forward problem for 2x2 mesh:
%--------------------------------------------------------------------------
% Author: Stefan Szyniszewski and Edward Street
% Date:   10/03/2021
%         14/06/2024
% Description: the script to solve the forward problem, i.e. computing
% displacement, strains, and stress for a problem with known boundary
% conditons, loads and elastic materials properties.
%
% Varargin:
%         'noplot'  = omits plotting the element displacements
%         'noprint' = omits printing to the console
%         'notime'  = omits computing the run time
%
% Example Usage:
%         [U,eps,sig] = forward_problem_2x2_elems_ES([25 25 25 25],[0.2 0.2 0.2 0.2],-3)
%
%--------------------------------------------------------------------------
% clearvars; close all; 
% clc;
% [thisPath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
% cd(thisPath);                                                             % change directory to current path
% addpath('functions');                                                       % add path to functions folder

%% Start clock
if ~any(strcmp(varargin,'notime'))
    tic;
end

%% Creating an element stiffness matrix for 2D bilinear elements in plane-stress conditions:
syms E nu x y real
% Plane-stress elastic stiffness matrix:
De = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
De_size = size(De); De_class = class(De);
% Shape functions of bilinear elements (assuming local and global
% coordinates are aligned, i.e. det([J])=0):
N1 = 1/4*(1-x)*(1-y); % bottom left
N2 = 1/4*(1+x)*(1-y); % bottom right
N3 = 1/4*(1+x)*(1+y); % top right
N4 = 1/4*(1-x)*(1+y); % top left
% Strain-displacement matrix:
B = [ diff(N1,x),          0, diff(N2,x),          0, diff(N3,x),          0, diff(N4,x),          0;
               0, diff(N1,y),          0, diff(N2,y),          0, diff(N3,y),          0, diff(N4,y);
      diff(N1,y), diff(N1,x), diff(N2,y), diff(N2,x), diff(N3,y), diff(N3,x), diff(N4,y), diff(N4,x) ];
% Element stiffness matrix:  
ke = int(int(B'*De*B,-1,1),-1,1);
ke_size = size(ke); ke_class = class(ke);
% Matrix after factorising common terms:
ke_fact = simplify(ke/E*24*(1-nu^2),'Steps',10);

%% Example construction of the global stiffness matrix for a 2x2 mesh:
%
%    ^ 2          ^ 8           ^14
%    |  1         |  7          |  13
%    o-->---------o-->----------o-->
% [1]|         [4]|          [7]|
%    |            |             |
%    |     (1)    |      (3)    |
%    ^ 4          ^ 10          ^ 16
%    |  3         |  9          |  15
%    o-->---------o-->----------o-->
% [2]|         [5]|          [8]|
%    |            |             |
%    |     (2)    |      (4)    |
%    ^ 6          ^ 12          ^ 18
%    |  5         |  11         |  17
%    o-->---------o-->----------o-->
% [3]          [6]           [9]
%
%  --> 1  - dof number
%  [1] o  - node
%  (1)    - elemenent
%

%% Generation of symbolic variables for each element:
E_g = sym('E_%d',[1 4],'real');                                                    % symbolic Young modulus in each element 1, 2, ..., 4
nu_g = sym('nu_%d',[1 4],'real');                                                  % symbolic Poisson ratio in each element 1, 2, ..., 4

%% Global stiffness matrix:
nelx = 2;                                                                   % number of elements in x-direction
nely = 2;                                                                   % number of elements in y-direction
neltot = nelx*nely;                                                         % total number of elements
ngnodes = (nelx+1)*(nely+1);                                                % total number of nodes
nldofs = 8;                                                                 % number of dof of each elements
ngdofs = 2*ngnodes;                                                         % total number of dofs (2 dof per each node)
gnodes = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);                        % Matrix with global node numbers                        
edofvec = reshape(2*gnodes(1:end-1,1:end-1)+1,neltot,1);
edofmat = repmat(edofvec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],neltot,1); % Matrix giving global dofs for each element (row)
%
%               [3	4	9	10	7	8	1	 2x ]  <- Element 1
%  edofMat =    [5	6	11	12	9	10	3	 4 ]  <- Element 2
%               [9	10	15	16	13	14	7	 8 ]  <- Element 3
%               [11	12	17	18	15	16	9	10 ]  <- Element 4
%
Kg = sym(zeros(ngdofs,ngdofs));
for i = 1:neltot                                                            % for each element:
    tmp = subs(ke,E,E_g(i));                                                % substitute E for E_g(i)
    ke_i = subs(tmp,nu,nu_g(i));                                            % substitute nu for nu_g(i)
    for j = 1:2:nldofs                                                      
        for k = 1:2:nldofs                                                  
            Kg(edofmat(i,j):edofmat(i,j+1),edofmat(i,k):edofmat(i,k+1)) = ...
                Kg(edofmat(i,j):edofmat(i,j+1),edofmat(i,k):edofmat(i,k+1)) + ...
                ke_i(j:j+1,k:k+1);
        end
    end
end
clear i j k tmp ke_i
Kg = simplify(Kg,'Steps',10);
Kg_size = size(Kg); Kg_class = class(Kg);

%% Load and boundary conditions (BC):
%   
%   || P          || P          || P
%   \/            \/            \/
%    ^ 2_         ^ 8_          ^14_
%    |  1         |  7          |  13
%    o-->---------o-->----------o-->
%    |            |             |
%    |            |             |
%    |     (1)    |      (3)    |
%    ^ 4          ^ 10          ^ 16
%    |  3         |  9          |  15
%    o-->---------o-->----------o-->
%    |            |             |
%    |            |             |
%    |     (2)    |      (4)    |
%    ^ 6_         ^ 12_         ^ 18_
%    |  5_        |  11         |  17
% |> o-->---------o-->----------o-->
%    /\           /\            /\
%    ""           ""            ""
force_dofs = [2 8 14];                                                      % dof with force application
fixed_dofs = [5 6 12 18];                                                   % dof with fixed dof 
all_dofs = 1:ngdofs;
free_dofs = setdiff(all_dofs,fixed_dofs);
%% Remove rows and columns corresponding to fixed BC
K_free = Kg(free_dofs,free_dofs);                                           % stiffness matrix corresponding to unknown displacements

%% Create loading vector:
Force = sym('F_%d',[ngdofs 1],'real');                                     % initialize symbolic force vector
syms P real                                                                 % symbolic load value
Force(free_dofs) = zeros(size(free_dofs));                                  % assign zeros to dofs without any load applied                                                    
Force(force_dofs) = P * ones(size(force_dofs,2),1);                         % assign loads to selected dofs
Force([force_dofs(1) force_dofs(end)]) = 0.5*P;
Force_size = size(Force); Force_class = class(Force);

%% Substituting specific values for variables (using numeric list):
E = E_arg(1);
nu = nu_arg(1);
P = P_arg; % -ve 'pushes into' the mesh
De_subs = subs(De);
ke_subs = subs(ke);                                                         % Symbolic expression using fractions
ke_subn = double(ke_subs);                                                  % Numerical expression in double precision
Force_subs = subs(Force);
E_num = E_arg;
nu_num = nu_arg;
for i = 1:4
    eval(sprintf('E_%d = E_num(%d);',i,i));
    eval(sprintf('nu_%d = nu_num(%d);',i,i));
end
Kg_subs = subs(Kg);                                                         % Symbolic expression using fractions
Kg_subn = double(Kg_subs);                                                  % Numerical expression in double precision

%% Solve for the global, unknown displacements:
% U(free_dofs) = K_free\Force(free_dofs); % Takes forever
% disp(U)
U = sym('U_%d',[ngdofs 1],'real');
for i = fixed_dofs
    eval(sprintf('U_%d = 0;',i))
end
U_subs = subs(U);
U_subs(free_dofs) = Kg_subs(free_dofs,free_dofs)\Force_subs(free_dofs);
U_subn = double(U_subs);
U_size = size(U); U_class = class(U);

%% Solve for strains in each element:
% use edofMat
for i = 1:4
    eval(sprintf('epsilon_%d = B*U_subs(edofmat(%d,:));',i,i))
end
epsilon_size = size(epsilon_1); epsilon_class = class(epsilon_1);

%% Solve for stresses in each element:
for i = 1:4
    eval(sprintf('sigma_%d = De_subs*epsilon_%d;',i,i))
end
sigma_size = size(sigma_1); sigma_class = class(sigma_1);

%% Print to console
if ~any(strcmp(varargin,'noprint'))
    fprintf('Elastic Stiffness Matrix for Plane Stress:\n');
    fprintf('========================\n');
    fprintf('De =\n');
    fprintf('%d x %d %s \n',De_size(1),De_size(2),De_class);
    disp(De);
    
    fprintf('Element stiffness matrix:\n');
    fprintf('========================\n');
    fprintf('ke =\n');
    fprintf('%d x %d %s \n',ke_size(1),ke_size(2),ke_class);
    disp(ke);
    
    fprintf('Global stiffness matrix:\n');
    fprintf('========================\n');
    fprintf('Kg =\n');
    fprintf('%d x %d %s \n',Kg_size(1),Kg_size(2),Kg_class);
    disp(Kg)
    
    fprintf('Load Vector:\n');
    fprintf('========================\n');
    fprintf('f =\n');
    fprintf('%d x %d %s \n',Force_size(1),Force_size(2),Force_class);
    disp(Force)
    
    fprintf('Substituting Variables:\n');
    fprintf('========================\n');
    fprintf('E = %f \t \x03BD = %f \t P = %f :\n\n',E,nu,P)
    fprintf('De =\n');
    fprintf('%d x %d %s \n',De_size(1),De_size(2),De_class);
    disp(De_subs);
    fprintf('ke =\n');
    fprintf('%d x %d %s \n',ke_size(1),ke_size(2),ke_class);
    disp(ke_subs);
    fprintf('Kg =\n');
    fprintf('%d x %d %s \n',Kg_size(1),Kg_size(2),Kg_class);
    disp(Kg_subs);
    
    fprintf('Solving Displacements:\n');
    fprintf('========================\n');
    fprintf('U =\n');
    fprintf('%d x %d %s \n',U_size(1),U_size(2),U_class);
    disp(U_subn);
    
    fprintf('Solving Strains:\n');
    fprintf('========================\n');
    for i = 1:4
        fprintf('epsilon_%d =\n',i);
        fprintf('%d x %d %s \n',epsilon_size(1),epsilon_size(2),epsilon_class);
        eval(sprintf('disp(epsilon_%d)',i));
    end
    
    fprintf('Solving Stresses:\n');
    fprintf('========================\n');
    for i = 1:4
        fprintf('sigma_%d =\n',i);
        fprintf('%d x %d %s \n',sigma_size(1),sigma_size(2),sigma_class);
        eval(sprintf('disp(sigma_%d)',i));
    end
end

%% Plotting the element displacements
if ~any(strcmp(varargin,'noplot'))
    coords_old = [0 4; 0 2; 0 0; 2 4; 2 2; 2 0; 4 4; 4 2; 4 0];
    coords_new = coords_old+reshape(U_subn,2,9)';
    stencil = [edofmat([1 2 3 4],2:2:end)/2 edofmat([1 2 3 4],2)/2];
    figure; axis equal; hold on;
    for i = 1:4
        for j = 1:4
            plot([coords_old(stencil(i,j),1) coords_old(stencil(i,j+1),1)],[coords_old(stencil(i,j),2) coords_old(stencil(i,j+1),2)],'-ok','LineWidth',1);
        end
    end
    for i = 1:4
        for j = 1:4
            plot([coords_new(stencil(i,j),1) coords_new(stencil(i,j+1),1)],[coords_new(stencil(i,j),2) coords_new(stencil(i,j+1),2)],'-or','LineWidth',1);
        end
    end
    legend([plot(nan,nan,'-ok'),plot(nan,nan,'-or')],["Before","After"])
end

%% Returning variables
U_return = U_subn;
Epsilon_return = [epsilon_1 epsilon_2 epsilon_3 epsilon_4];
Sigma_return = [sigma_1 sigma_2 sigma_3 sigma_4];
% Integrate stress and strain across each element
Epsilon_return = double(int(int(Epsilon_return,-1,1),-1,1));
Sigma_return = double(int(int(Sigma_return,-1,1),-1,1));

%% Stop clock
if ~any(strcmp(varargin,'notime'))
    toc;
end

end