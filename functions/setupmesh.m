function [mesh] = setupmesh(lx,ly,nels_x,nels_y,pressure,E_mat,v_mat,bc_type)

%Mesh generation and input information
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   07/11/2019
% Description:
% Finite element setup information for a rectangular domain discretised
% with 4-noded elements. 
%
%--------------------------------------------------------------------------
% [coord,etpl,fext,bc,ngp,E,v] = SETUPMESH(lx,ly,nels_x,nels_y,pressure,E_init,v_init)
%--------------------------------------------------------------------------
% Input(s):
%           - E_init - Young modulus for the initial mesh
%           - v_init - Poisson ratio for the initial mesh
%           - lx    - domain size in x-direction
%           - ly    - domain size in y-direction
%           - nels_x - number of elements in x-direction
%           - nels_y - number of elements in y-direction
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

% Boundary condition type                       Flag (BCt)
% 
% homogeneous Dirichlet (u = 0)                     1
% inhomogeneous Dirichlet (u \neq 0 )              -1
% homogeneous Neumann (p = 0)                       2
% inhomogeneous Neumann (p \neq 0)                 -2
% mixed (roller)                                    3
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

%% Boundary condition information 
if strcmp(bc_type,'roller')
    BCt    = [2 -2  2 3];                                                   % boundary condition flags
    BCp    = [0 -pressure 0 0];                                             % boundary condition pressures (normal) for each edge
                                                                            % Second entry is pressure applied to face (2)
elseif strcmp(bc_type,'fixed')
    BCt    = [2 -2  2 1];                                                          
    BCp    = [0 -pressure 0 0];                                             % boundary condition pressures (normal) for each edge
elseif strcmp(bc_type,'pressure')                                           % Negative value for compression and positive for tension
%     BCt    = [2 -2  3 -2]; 
    BCt    = [2 -2  2 -2]; 
    BCp    = [0 -pressure 0 -pressure];                                     % boundary condition pressures (normal) for each edge

else
    BCt    = [2 -2  2 1];  
    BCp    = [0 -pressure 0 0];  
end
    
BCu    = [0 0;                                                              % boundary condition displacements (x,y) for each edge
          0 0; 
          0 0; 
          0 0];

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

[etpl,coord] = formMesh(nels_x,nels_y,lx,ly);                               % element topology and nodal coordinates

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
        
        p = BCp(edge)*lx/nels_x;                                            % load over each element segment 
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

tol    = 1e-9;                                                              % numerical tolerance
if strcmp(bc_type,'pressure')  
%
%               ^ y
%               |
%               |
%                           (2)
%                ---------------------------
%               |            o|             |
%               |                           |
%               |                           |
%               |                           |
%           (1) |            o              | (3)
%               |            -              |
%               |                           |
%               |                           |
%               |                           |
%               |            o|             |
%                ---------------------------      -----> x
%                           (4)
%
%
% The problem will only be symmetric if an even number of elements are 
% used in the x and y directions.  If not, the boundary conditions will be 
% imposed off centre. 

    % Roller boundary conditions applied at x position lx/2.
    % Fixed boundary condition applied at (x,y) = (lx/2,ly/2).
    % These are implemented below.
    hx = lx/nels_x;                                                         % element spacing (x direction)
    hy = ly/nels_y;                                                         % element spacing (y direction)
    for node = 1:nodes                                                      % loop over all nodes
        if abs(coord(node,1)-lx/2)<tol || abs(coord(node,1)-lx/2-hx/2)<tol  % if x position located at lx/2
%             bc(node*nD-1,:) = [node*nD-1 0];                                % fix x direction
            if abs(coord(node,2)-ly/2)<tol || abs(coord(node,2)-ly/2-hy/2)<tol  % if y position located at ly/2
               % bc(node*nD-1,:) = [node*nD-1 0];                            % fix x direction
                bc(node*nD,:) = [node*nD   0];                              % fix y direction 
            end

            if abs(coord(node,2)-0)<tol || abs(coord(node,2)-0-hy/2)<tol    % if y position located at 0
                bc(node*nD-1,:) = [node*nD-1 0];                            % fix x direction
            end

            if abs(coord(node,2)-ly)<tol || abs(coord(node,2)-ly-hy/2)<tol    % if y position located at ly
                bc(node*nD-1,:) = [node*nD-1 0];                            % fix x direction
            end

        end
    end
end


bc = bc(bc(:,1)>0,:);                                                       % remove unused bc rows

if (isscalar(E_mat))
    E = E_mat*ones(nels,1);                                                 % propagate scalar into all cells of the vector  
else 
    E = E_mat;                                                              % use vector input
end

if (isscalar(v_mat))
    v = v_mat*ones(nels,1);                                                 % propagate scalar into all cells
else
    v = v_mat;                                                              % use vector input
end
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
