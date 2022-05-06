% PURPOSE
%     Symbolic solution for 2x2 elements under load, with known BC
%
% DEPENDENCIES:
% 
% RELATED SCRIPTS:
%     stiffness_matrices.m
% Date:
%     23-March-2022
%  ----------------------------------------------------------------
clearvars; close all; clc;
[thisPath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(thisPath); %addpath('functions') 

%% Force-displacement system of equation
% Formulate element stiffness matrix first, and assemble the global
% stiffness matrix:

%% Element stiffness matrix for 2D bilinear elements in plane-stress conditions:
syms E nu x y                                                               % introducte symbolic variables
% Plane-stress elastic stress-strain (material stiffness) matrix:
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
% Element stiffness matrix  
ke = int(int(B'*De*B,-1,1),-1,1);
disp(ke)
% Matrix after factorising common terms
ke_fact = simplify(ke/E*24*(1-nu^2));
disp( E/(24*(1-nu^2)) )
disp(ke_fact)

%% Construction of the global stiffness matrix for a 2x2 mesh:
% ASCII art?

nelx = 2;                                                                   % number of elements in x-direction
nely = 2;                       
neltot = nelx*nely;                                                         % total number of elements
ngnodes = (nelx+1)*(nely+1);                                                
nldofs = 8;                                                                 % 
ngdofs = 2*ngnodes;
gnodes = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofvec = reshape(2*gnodes(1:end-1,1:end-1)+1,neltot,1);
edofmat = repmat(edofvec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],neltot,1);
K = sym(zeros(ngdofs,ngdofs));
for i = 1:neltot
    for j = 1:2:nldofs
        for k = 1:2:nldofs
            K(edofmat(i,j):edofmat(i,j+1),edofmat(i,k):edofmat(i,k+1)) = K(edofmat(i,j):edofmat(i,j+1),edofmat(i,k):edofmat(i,k+1)) + ke(j:j+1,k:k+1);
        end
    end
end
disp(K)





%% Solve for the global, unknown displacements:


%% Solve for strains in each element:
% use edofMat


%% Solve for stresses in each element:
