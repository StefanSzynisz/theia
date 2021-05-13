 
function [uvw] = linearSolver(K,bc,fext)

%Linear solver 
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   07/11/2019
% Description:
% Function to determine the displacements given the stiffness matrix,
% external forces and boundary conditions.
%
%--------------------------------------------------------------------------
% [uvw] = LINEARSOLVER(K,bc,fext,nDoF)
%--------------------------------------------------------------------------
% Input(s):
% K     - global stiffness matrix (sparse)
% bc    - matrix of boundary conditions [degree of freem, displacement]
% fext  - external force vector
%--------------------------------------------------------------------------
% Ouput(s);
% uvw   - global displacement vector
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------

nDoF = length(fext);                                                        % number of degrees of freedom
fd=(1:nDoF);                                                                % free degrees of freedom
uvw=zeros(nDoF,1);                                                          % zero displacements
if length(bc)>1
    fd(bc(:,1))=[];                                                         % remove constrained degrees of freedom
    uvw(bc(:,1))=bc(:,2);                                                   % boundary conditions
    fext(fd)=fext(fd)-K(fd,bc(:,1))*uvw(bc(:,1));                           % modify external forces for imposed displacements
end
uvw(fd) = K(fd,fd)\fext(fd);                                                % linear solve
fprintf('%s\n','   system solved');