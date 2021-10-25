function [sig,strain,uv] = LEfe(mesh)

%LEfe - Linear Elastic finite element solver
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   02/04/2020
% Description:
% Finite element solver for 2D plane stress analysis based on a vectorised
% implementation. 
%
%--------------------------------------------------------------------------
% [sig,strain] = LEFE(mesh)
%--------------------------------------------------------------------------
% Input(s):
% mesh   - structured array of mesh data
% itnum  - iteration number 
%--------------------------------------------------------------------------
% Ouput(s);
% sig     - Cauchy stress at the centre of each element (nels,3)
% strain  - strain at the centre of each element (nels,3)
% 
%--------------------------------------------------------------------------
% See also:
% 
% FORMKVEC            - global stiffness formation (vectorised)
% AVBCSVEC            - average boundary conditions (vectorised)
% LINEARSOLVER        - linear solver
% POSTPRO             - Cauchy stress determination
% MAKEVTK             - VTK mesh file generation
% MAKEVTKMP           - VTK stress data file generation
%--------------------------------------------------------------------------
K    = formKvec(mesh);                                                      % vectorised stiffness calculation
fext = mesh.fext;                                                           % external force vector

uvw = linearSolver(K,mesh.bc,fext);                                         % displacement solution
[sig,strain,uv,xy] = postPro(mesh,uvw);                                     % vectorised stress/displacement calculation

