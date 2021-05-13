function [K] = formKvec(mesh)

%Global stiffness matrix calculation
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   01/04/2020
% Description:
% Function to determine the global stiffness matrix based on a vectorised
% algorithm.  The algorithm allows for different stiffness values in
% different elements. 
%
%--------------------------------------------------------------------------
% [K] = FORMKVEC(coord,etpl,ngp,E,v,nelblo)
%--------------------------------------------------------------------------
% Input(s):
% mesh   - structured array of finite element data.  The function requires
%          the following:
%           - coord  - nodal coordinates (nodes,nD)
%           - etpl   - element topology (nels,nen)
%           - ngp    - number of Gauss points (scalar)
%           - E      - Young's modulus of each element (nels,1)
%           - v      - Poisson's ratio of each element (nels,1)
%           - nelblo - number of elements per block (scalar)
%--------------------------------------------------------------------------
% Ouput(s);
% K      - global stiffness matrix (sparse)
%--------------------------------------------------------------------------
% See also:
%
% DERSF2D   - shape function derivatives & Gauss weights
%--------------------------------------------------------------------------

coord  = mesh.coord;                                                        % nodal coordinates
etpl   = mesh.etpl;                                                         % element topology
ngp    = mesh.ngp;                                                          % number of Gauss points
E      = mesh.E;                                                            % Young's modulus
v      = mesh.v;                                                            % Poisson's ratio
nelblo = mesh.nelblo;                                                       % number of elements per block

[nels,nen]=size(etpl);                                                      % no. elements and nodes/element
[~,nD]=size(coord);                                                         % no. nodes and dimensions
nedof=nen*nD;                                                               % no. element degrees of freedom
[wp,dNr]=derSF2D(ngp,nen);                                                  % Gauss weights and shape function derivatives

C0=E./(1-v.^2); C11=C0*1; C12=C0.*v; C33=C0.*(1-v)/2;                       % elastic constants (plane stress)

nelblo=min(nelblo,nels);                                                    % check on block size < no. elements
nbloc=ceil(nels/nelblo);                                                    % no. blocks
Kall=zeros(nedof*(nedof+1)/2,nels);                                         % global stiffness entries
Kb=zeros(nelblo,nedof*(nedof+1)/2);                                         % block stiffness entries
invJx=zeros(nelblo,nD);                                                     % Jacobian inverse entries (x-direction)
invJy=zeros(nelblo,nD);                                                     % Jacobian inverse entries (y-direction)

ilow=1; iup=nelblo;                                                         % block range for first loop
for ib=1:nbloc                                                              % block loop 
    coordx=reshape(coord(etpl(ilow:iup,:),1), nelblo, nen);                 % x coordinates
    coordy=reshape(coord(etpl(ilow:iup,:),2), nelblo, nen);                 % y coordinates
    
    C1 = C11(ilow:iup);                                                     % elastic properties xx-xx, yy-yy
    C2 = C12(ilow:iup);                                                     % elastic properties xx-yy, yy-xx
    C3 = C33(ilow:iup);                                                                              % elastic properties xy-xy
    
    Kb(:)=0;                                                                % zero block stiffness 
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
        x=1;        
        for i = 1:nen                                                       % stiffness matrix entries
          for j = i:nen
            Kb(:,x) = Kb(:,x)+(C1.*dNx(:,i).*dNx(:,j)+C3.*dNy(:,i).*dNy(:,j)).*w; x=x+1;
            Kb(:,x) = Kb(:,x)+(C2.*dNx(:,i).*dNy(:,j)+C3.*dNy(:,i).*dNx(:,j)).*w; x=x+1;
          end
          
          for j = i:nen
            if(j>i)
              Kb(:,x)=Kb(:,x)+(C2.*dNy(:,i).*dNx(:,j)+C3.*dNx(:,i).*dNy(:,j)).*w; x=x+1;
            end
            Kb(:,x) = Kb(:,x)+(C1.*dNy(:,i).*dNy(:,j)+C3.*dNx(:,i).*dNx(:,j)).*w; x=x+1;
          end

        end
        Kall(:,ilow:iup)=Kb';                                               % store block stiffness enrties
    end
    
    ilow=ilow+nelblo;                                                       % block index adjustment
    if(ib==nbloc-1)                                                         % final block adjustment
        nelblo=nels-iup;
        Kb=zeros(nelblo, nedof*(nedof+1)/2);
        invJx=zeros(nelblo, nD);
        invJy=zeros(nelblo, nD);
    end
    iup=iup+nelblo;                                                         % update block range
end

elDoF = zeros(nedof, nels);                                                 % element degrees of freedom
elDoF(1:nD:end,:)=(nD*(etpl-1)+1)';
elDoF(2:nD:end,:)=(nD*(etpl-1)+2)';
indx_j=repmat(int32(1:nedof),nedof,1); 
indx_i=indx_j';
indx_i=tril(indx_i); indx_i=indx_i(:); indx_i=indx_i(indx_i>0);
indx_j=tril(indx_j); indx_j=indx_j(:); indx_j=indx_j(indx_j>0);
K_i=elDoF(indx_i(:),:);                                                     % row entry locations
K_j=elDoF(indx_j(:),:);                                                     % column entry locations
indx=K_i < K_j;
tmp=K_j(indx);
K_j(indx)=K_i(indx);
K_i(indx)=tmp;
K=sparse(K_j,K_i,Kall);                                                     % sparse assembly
K=triu(K)+triu(K,1)';                                                       % full matrix assembly
fprintf('%s\n','   stiffness assembled');