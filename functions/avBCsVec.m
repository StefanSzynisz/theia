function [K] = avBCsVec(mesh)

%Average BC stiffness calculation (vectorised)
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   04/04/2020
% Description:
% Function to determine the global stiffness matrix additions for average
% bounadry conditions for two-dimensional analysis. 
%
%--------------------------------------------------------------------------
% [K] = AVBCSVEC(coord,etpl,ngp,E,v,nelblo)
%--------------------------------------------------------------------------
% Input(s):
% mesh   - structured array of finite element data.  The function requires
%          the following:
%           - coord  - nodal coordinates (nodes,nD)
%           - etpl   - element topology (nels,nen)
%           - ngp    - number of Gauss points (scalar)
%           - nelblo - number of elements per block (scalar)
%--------------------------------------------------------------------------
% Ouput(s);
% K      - global stiffness matrix additional components (nDoF,3)
%--------------------------------------------------------------------------
% See also:
%
% DERSF2D       - shape function derivatives & Gauss weights
% SF            - shape functions
%--------------------------------------------------------------------------

coord  = mesh.coord;                                                        % nodal coordinates
etpl   = mesh.etpl;                                                         % element topology
ngp    = mesh.ngp;                                                          % number of Gauss points
nelblo = mesh.nelblo;                                                       % number of elements per block

[nels,nen] = size(etpl);                                                    % no. elements and nodes/element
[nodes,nD] = size(coord);                                                   % no. nodes and dimensions
[wp,dNr]   = derSF2D(ngp,nen);                                              % Gauss weights and shape function derivatives
[N]        = SF(1:ngp,ngp,nen);                                             % Shape functions 
nelblo=min(nelblo,nels);                                                    % check on block size < no. elements
nbloc=ceil(nels/nelblo);                                                    % no. blocks
invJx=zeros(nelblo,nD);                                                     % Jacobian inverse entries (x-direction)
invJy=zeros(nelblo,nD);                                                     % Jacobian inverse entries (y-direction)
Kr = zeros(nodes*nD,1);                                                     % zero stiffness storage
Ku = zeros(nodes*nD,1);
Kv = zeros(nodes*nD,1);
ilow=1; iup=nelblo;                                                         % block range for first loop
for ib=1:nbloc                                                              % block loop 
    coordx=reshape(coord(etpl(ilow:iup,:),1), nelblo, nen);                 % x coordinates
    coordy=reshape(coord(etpl(ilow:iup,:),2), nelblo, nen);                 % y coordinates
    for gp=1:ngp                                                            % Gauss point loop
        gpi=int32(nD*gp-(nD-1:-1:0));                                       % index for Gauss point entries
        dNdui=dNr(gpi,:);                                                   % local derivatives
        Jx=coordx*dNdui.';                                                  % Jacobian entries                                  
        Jy=coordy*dNdui.';
        detJ=Jx(:,1).*Jy(:,2)-Jx(:,2).*Jy(:,1);                             % Jacobian determinant
        invdetJ=1./detJ;                                                    % inverse Jac. det.                                                  
        invJx(:,1)= Jy(:,2).*invdetJ;                                       % inverse Jac.
        invJx(:,2)=-Jy(:,1).*invdetJ;
        invJy(:,1)=-Jx(:,2).*invdetJ;
        invJy(:,2)= Jx(:,1).*invdetJ;
        dNx=invJx*dNdui;                                                    % global shape function derivatives (x) 
        dNy=invJy*dNdui;                                                    % global shape function derivatives (y) 
        w=detJ*wp(gp);                                                      % weight inc. det. Jac.
        for i = 1:nen            
            eDoFs = nD*(etpl(ilow:iup,i)-1);                                % degrees of freedom of stiffness entries
            Ku(eDoFs+1) = Ku(eDoFs+1) + N(i,gp).*w;                         % x displacement restraint
            Kv(eDoFs+2) = Kv(eDoFs+2) + N(i,gp).*w;                         % y displacement restraint  
            Kr(eDoFs+1) = Kr(eDoFs+1) - dNy(:,i).*w;                        % rotation restraint (x direction)
            Kr(eDoFs+2) = Kr(eDoFs+2) + dNx(:,i).*w;                        % rotation restraint (y direction)
        end
    end
    ilow=ilow+nelblo;                                                       % block index adjustment
    if(ib==nbloc-1)                                                         % final block adjustment
        nelblo=nels-iup;
        invJx=zeros(nelblo, nD);
        invJy=zeros(nelblo, nD);
    end
    iup=iup+nelblo;                                                         % update block range
end
K = [Ku Kv Kr];
fprintf('%s\n','   average boundary conditions assembled');