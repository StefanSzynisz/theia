% forward problem for 2x2 mesh:
%--------------------------------------------------------------------------
% Author: Stefan Szyniszewski and Edward Street
% Date:   10/03/2021
% Description: the script to solve the forward problem, i.e. computing
% displacement, strains, and stress for a problem with known boundary
% conditons, loads and elastic materials properties.
%
%--------------------------------------------------------------------------
clearvars; close all; clc;
% [thisPath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
% cd(thisPath);                                                             % change directory to current path
addpath('functions');                                                       % add path to functions folder
clear; tic; 

%% Constructing the element stiffness matrix

% Creating an element stiffness matrix for 2D bilinear elements in plane-stress conditions
syms E nu x y
% Plane-stress elastic stiffness matrix
De = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
% Shape functions of bilinear elements (assuming local and global coordinates are aligned)
N1 = 1/4*(1-x)*(1-y);
N2 = 1/4*(1+x)*(1-y);
N3 = 1/4*(1+x)*(1+y);
N4 = 1/4*(1-x)*(1+y);
% Strain-displacement matrix
B = [ diff(N1,x),          0, diff(N2,x),          0, diff(N3,x),          0, diff(N4,x),          0;
               0, diff(N1,y),          0, diff(N2,y),          0, diff(N3,y),          0, diff(N4,y);
      diff(N1,y), diff(N1,x), diff(N2,y), diff(N2,x), diff(N3,y), diff(N3,x), diff(N4,y), diff(N4,x) ];
% Element stiffness matrix  
ke = int(int(B'*De*B,-1,1),-1,1);
disp(ke)
% Matrix after factorising common terms
ke_fact = simplify(ke/E*24*(1-nu^2),'Steps',10);

%% Evaluating matrix slices
A = sym(zeros(size(ke_fact)));
B = sym(zeros(size(ke_fact)));
dims = size(ke_fact);
for i = 1:dims(1)
    for j = 1:dims(2)
        coeffs = sym2poly(ke_fact(i,j));
        B(i,j) = coeffs(1);
        A(i,j) = coeffs(2);
    end
end
A11 = A(1:dims(1)/2,1:dims(2)/2);
A12 = A(1:dims(1)/2,dims(2)/2+1:dims(2));
B11 = B(1:dims(1)/2,1:dims(2)/2);
B12 = B(1:dims(1)/2,dims(2)/2+1:dims(2));
disp(A11)
disp(A12)
disp(B11)
disp(B12)
%% Substituting specific values for variables
E = 1;
nu = 0.3;
ke_subs = subs(ke); % Symbolic expression using fractions
ke_subn = double(ke_subs); % Numerical expression in double precision
disp(ke_subs)
disp(ke_subn)
% Example construction of the global stiffness matrix for a 2x2 mesh
nelx = 2;
nely = 2;
neltot = nelx*nely;
ngnodes = (nelx+1)*(nely+1);
nldofs = 8;
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
