function [mesh] = setupmesh

%Mesh generation and input information
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   07/11/2019
% Description:
% Finite element setup information for a rectangular domain discretised
% with 4-noded elements. 
%
% Boundary condition type                       Flag (BCt)
% 
% homogeneous Dirichlet (u = 0)                     1
% inhomogeneous Dirichlet (u \neq 0                -1
% homogeneous Neumann (p = 0)                       2
% inhomogeneous Neumann (p \neq 0)                 -2
% mixed (roller)                                    3
%
%
% Boundary ordering
%
%               ^ y
%               |
%               |
%                           (2)
%                -------------------------
%               |                         |
%               |                         |
%               |                         |
%               |                         |
%           (1) |                         | (3)
%               |                         |
%               |                         |
%               |                         |
%               |                         |
%               |                         |
%                -------------------------      -----> x
%                           (4)
%
% Displacement directions are defined in terms of the global coordinate
% directions.
% 
% Pressures are relative to the outward normal direction to the boundary 
% (tensile pressures taken as positive).
%
% If no Dirichlet (displacement) boundary conditions are applied the
% average displacement boundary flag, avBC, should be set to 1. Note that
% the average bounadry condition should only be included if required,
% including the condition when other displacement boundary conditions are
% imposed will result in spurious results as the displacement of the nodes
% will be over constrained.
%
%--------------------------------------------------------------------------
% [coord,etpl,fext,bc,ngp,E,v] = SETUPMESH
%--------------------------------------------------------------------------
% Input(s):
% 
%--------------------------------------------------------------------------
% Ouput(s);
% mesh  - structured array with finite element data, including:
%           - coord  - nodal coordinates
%           - etpl   - element topology
%           - fext   - external force vector
%           - bc     - displacement boundary conditions
%           - ngp    - number of Gauss points (total per element)
%           - E      - Young's modulus 
%           - v      - Poissons ratio
%           - ngpP   - number of Gauss points for plotting
%           - VTKout - VTK output flag 
%           - nelblo - no. elements/block (for vectorisation)
%--------------------------------------------------------------------------
% See also:
% 
% FORMMESH   - 2D mesh generation
% NODES2DOFS - nodes to degrees of freedom
%--------------------------------------------------------------------------

%% Basic material properties and domain size
E      = 200e9;                                                             % Young's modulus     
v      = 0.3;                                                               % Poissons ratio

lx     = 1;                                                                 % domain length in the x direction
ly     = 1;                                                                 % domain length in the y direction

%% Numbers of elements
nelsx = 20;                                                                 % number of elements in the x direction
nelsy = nelsx;                                                              % number of elements in the y direction

%% Boundary condition information 
BCt    = [2 -2  2 3];                                                     % boundary condition flags
BCu    = [0 0;                                                              % boundary condition displacements (x,y) for each edge
          0 0; 
          0 0; 
          0 0];                                                             
BCp    = [0 -1e9 0 0];                                           % boundary condition pressures (normal) for each edge

%% Other settings/output options
VTKout = 1;                                                                 % VTK output flag 
nelblo = 512;                                                               % no. elements/block (for vectorisation)

%--------------------------------------------------------------------------   do not change below this line

ngp    = 4;                                                                 % number of Gauss points per element
ngpP   = 1;                                                                 % number of Gauss points for plotting (set to 1)

noE = 4;                                                                    % number of edges
BCnorm = [-1 0; 0 1; 1 0; 0 -1];                                            % boundary condition normal directions
BCxy   = [2 1 2 1];                                                         % boundary condition coordinate limit direction
BClims = [0  0  0  ly;                                                      % boundary limits
          0  lx ly ly;
          lx lx 0  ly;
          0  lx 0  0]; 

[etpl,coord] = formMesh(nelsx,nelsy,lx,ly);                                 % element topology and nodal coordinates

[nodes,nD] = size(coord);                                                   % number of nodes and dimensions
[nels,~] = size(etpl);                                                      % number of elements

                                                                          
bc   = zeros(nodes*nD,2);                                                   % boundary condition matrix
fext = zeros(nodes*nD,1);                                                   % external force vector
nodeList = (1:nodes)';                                                      % list of nodes in the analysis

for edge = 1:noE                                                            % loop over mesh edges
    
    BCe = BCt(edge);                                                        % BC edge type
    
    x = BClims(edge,1:2);                                                   % edge x limits
    y = BClims(edge,3:4);                                                   % edge y limits
    
    nn = coord(:,1)>=x(1) & coord(:,1)<=x(2) &... 
         coord(:,2)>=y(1) & coord(:,2)<=y(2);                               % logical flags of nodes on edge
    nn = nodeList(nn);                                                      % nodes on edge
    
    dofs  = nodes2dofs(nn,nD)';                                             % degrees of freedom for nodes on the edge
    dofsX = dofs(1:nD:end-1);                                               % x degrees of freedom
    dofsY = dofs(2:nD:end  );                                               % y degrees of freedom
    
    if abs(BCe) == 1                                                        % displacement (Dirichlet) boundary condition
        
        bc(dofsX,:) = [dofsX BCu(edge,1)*ones(size(dofsX))];                % x direction boundary conditions
        bc(dofsY,:) = [dofsY BCu(edge,2)*ones(size(dofsY))];                % y direction boundary conditions
        
    elseif BCe == -2                                                        % traction (Neumann) boundary conditions
        
        p = BCp(edge)*lx/nelsx;                                             % load over each element segment 
        p = p*ones(length(nn),1)* BCnorm(edge,:);                           % nodal forces [x y]
        p = reshape(p',length(dofs),1);                                     % nodal forces (vector format)
        if BCxy(edge) == 1
            i = coord(nn,1)==x(1) | coord(nn,1)==x(2);                      % logical flag for nodes at end of edge
        else 
            i = coord(nn,2)==y(1) | coord(nn,2)==y(2);                      % logical flag for nodes at end of edge
        end
        ii = (1:length(nn))';                                               % list [1:nn]
        i = nodes2dofs(ii(i),nD);                                           % corner dofs
        p(i) = p(i)/2;                                                      % reduce corner loads
        fext(dofs) = fext(dofs) + p;                                        % external forces
    
    elseif BCe == 3                                                         % mixed (roller) boundary condition 
        
        if BCxy(edge) == 1                                                  % x edge
            bc(dofsY,:) = [dofsY BCu(edge,2)*ones(size(dofsY))];            % constain y direction
            bc(dofsX(11),:) = [dofsX(11) BCu(edge,1)*ones(size(dofsX(11)))];
        else                                                                % y edge
            bc(dofsX,:) = [dofsX BCu(edge,1)*ones(size(dofsX))];            % constrain x direction
        end
        
    end
end
bc = bc(bc(:,1)>0,:);                                                       % remove unused bc rows

E = E*ones(nels,1);                                                         % vector of Young's modulus                       	
v = v*ones(nels,1);                                                         % vector of Poisson's ratios
%% Mesh data structure generation
mesh.etpl   = etpl;                                                          % element topology
mesh.coord  = coord;                                                         % nodal coordinates
mesh.fext   = fext;                                                          % external forces
mesh.bc     = bc;                                                            % boundary conditions
mesh.E      = E;                                                             % Young's modulus
mesh.v      = v;                                                             % Poisson's ratio
mesh.ngp    = ngp;                                                           % number of Gauss points
mesh.ngpP   = ngpP;                                                          % number of Gauss points for plotting
mesh.VTKout = VTKout;                                                        % VTK output flag 
mesh.nelblo = nelblo;                                                        % no. elements/block (for vectorisation)


function [dofs] = nodes2dofs(nn,nD)

%Node numbers to degrees of freedom
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   07/11/2019
% Description:
% Function to convert a list of nodes to a list of degrees of freedom. 
%
%--------------------------------------------------------------------------
% [dofs] = NODES2DOFS(nn,nD)
%--------------------------------------------------------------------------
% Input(s):
% nn     - list of nodes
% nD     - number of dimensions
%--------------------------------------------------------------------------
% Ouput(s):
% dofs   - degrees of freedom list
%--------------------------------------------------------------------------
% See also:
% 
%--------------------------------------------------------------------------

n = length(nn);                                                             % no. nodes
dofs=reshape(ones(nD,1)*nn'*nD-(nD-1:-1:0).'*ones(1,n),1,n*nD);             % degrees of freedom
