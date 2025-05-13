% Forward problem for 2x2 mesh:
%--------------------------------------------------------------------------
% Author: Stefan Szyniszewski and Edward Street
% Date:   10/03/2021
% Description: the script to solve the forward problem, i.e. computing
% displacement, strains, and stress for a problem with known boundary
% conditons, loads and elastic materials properties.
%
%--------------------------------------------------------------------------
clearvars; close all; clc;
[thisPath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(thisPath);                                                               % change directory to current path
addpath('functions');                                                       % add path to functions folder
clear; tic; 

%% Creating an element stiffness matrix for 2D bilinear elements in plane-stress conditions:
syms E nu x y
% Plane-stress elastic stiffness matrix:
De = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
% Shape functions of bilinear elements (assuming local and global
% coordinates are aligned):
N1 = 1/4*(1-x)*(1-y);
N2 = 1/4*(1+x)*(1-y);
N3 = 1/4*(1+x)*(1+y);
N4 = 1/4*(1-x)*(1+y);
% Strain-displacement matrix:
B = [ diff(N1,x),          0, diff(N2,x),          0, diff(N3,x),          0, diff(N4,x),          0;
               0, diff(N1,y),          0, diff(N2,y),          0, diff(N3,y),          0, diff(N4,y);
      diff(N1,y), diff(N1,x), diff(N2,y), diff(N2,x), diff(N3,y), diff(N3,x), diff(N4,y), diff(N4,x) ];
% Element stiffness matrix:  
disp('Element stiffness matrix:');
fprintf('======================== \n');
ke = int(int(B'*De*B,-1,1),-1,1);
disp(ke)
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
E_g = sym('E_%d',[1 4]);                                                    % symbolic Young modulus in each element 1, 2, ..., 4
nu_g = sym('nu_%d',[1 4]);                                                  % symbolic Poisson ratio in each element 1, 2, ..., 4

%% Globa stiffness matrix:
disp('Global stiffness matrix:');
fprintf('======================== \n');
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
%               [3	4	9	10	7	8	1	 2 ]  <- Element 1
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
disp(Kg)

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
all_dofs = [1:1:ngdofs];
free_dofs = setdiff(all_dofs,fixed_dofs);
%% Remove rows and columns corresponding to fixed BC
K_free = Kg(free_dofs,free_dofs);                                           % stiffness matrix corresponding to unknown displacements

%% Create loading vector:
disp('Global stiffness matrix:');
fprintf('======================== \n');
Force = sym('F_%d',[1 ngdofs])';                                            % initialize symbolic force vector
syms P                                                                      % symbolic load value
Force(free_dofs) = zeros(size(free_dofs));                                  % assign zeros to dofs without any load applied                                                    
Force(force_dofs) = P * ones(size(force_dofs,2),1);                         % assign loads to selected dofs
disp(Force)
%% Solve for the global, unknown displacements:



%% Substituting specific values for variables (using numeric list):
E_num = [25 25 25 25];
nu_num = [0.2 0.2 0.2 0.2];
% E = 1;
% nu = 0.3;
% ke_subs = subs(ke);                                                         % Symbolic expression using fractions
% ke_subn = double(ke_subs);                                                  % Numerical expression in double precision
% disp(ke_subs)
% disp(ke_subn)

% u_num;

%% Solve for strains in each element:
% use edofMat


%% Solve for stresses in each element:



