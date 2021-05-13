function [sig,strain,uv,xy] = postPro(mesh,uvw)

%Stress at the centre of each element 
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   07/11/2019
% Description:
% Function to determine the stress at the centre of each element using a
% vectorised calculation method.
%
%--------------------------------------------------------------------------
% [sig] = POSTPRO(coord,etpl,ngp,E,v,uvw,nelblo)
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
% uvw    - global displacement vector
%--------------------------------------------------------------------------
% Ouput(s);
% sig    - Cauchy stress in each element [sigXX sigYY sigXY]
% strain - strain in each element [epsXX epsYY espXY]
% uv     - displacement at the centre of each element [u v]
% xy     - centre location of each element [x y]
%--------------------------------------------------------------------------
% See also:
%
% DERSF2D   - shape function derivatives & Gauss weights
%--------------------------------------------------------------------------

coord  = mesh.coord;                                                        % nodal coordinates
etpl   = mesh.etpl;                                                         % element topology
ngp    = mesh.ngpP;                                                         % number of Gauss points
E      = mesh.E;                                                            % Young's modulus
v      = mesh.v;                                                            % Poisson's ratio
nelblo = mesh.nelblo;                                                       % number of elements per block

[nels,nen]=size(etpl);                                                      % no. elements and nodes/element
[~,nD]=size(coord);                                                         % no. nodes and dimensions 

[~,dNr]=derSF2D(ngp,nen);                                                   % Gauss weights and shape function derivatives

C0=E./(1-v.^2); C11=C0*1; C12=C0.*v; C33=C0.*(1-v)/2;                       % elastic constants

nelblo=min(nelblo,nels);                                                    % check on block size < no. elements
nbloc=ceil(nels/nelblo);                                                    % no. blocks
invJx=zeros(nelblo,nD);                                                     % Jacobian inverse entries
invJy=zeros(nelblo,nD);

sig   = zeros(nels,3);                                                      % zero stress storage
strain = zeros(nels,3);                                                     % zero stress storage

ilow=1; iup=nelblo;
for ib=1:nbloc                                                              % block loop 
    coordx=reshape(coord(etpl(ilow:iup,:),1), nelblo, nen);                 % x coordinates
    coordy=reshape(coord(etpl(ilow:iup,:),2), nelblo, nen);                 % y coordinates
    
    u = reshape(uvw(etpl(ilow:iup,:)*nD-1), nelblo, nen);                   % x displacement
    v = reshape(uvw(etpl(ilow:iup,:)*nD-0), nelblo, nen);                   % y displacement
    
    C1 = C11(ilow:iup);                                                     % block stiffness terms
    C2 = C12(ilow:iup);
    C3 = C33(ilow:iup);
       
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
        
        dNx=invJx*dNdui;                                                    % global shape function derivatives
        dNy=invJy*dNdui;
        
        epsXX = sum(dNx.*u,2);                                              % strains
        epsYY = sum(dNy.*v,2);
        epsXY = sum(dNx.*v+dNy.*u,2);
        
        sigXX = C1.*epsXX + C2.*epsYY;                                      % stresses
        sigYY = C2.*epsXX + C1.*epsYY;
        sigXY = C3.*epsXY;
        
    end
    
    sig(ilow:iup,1) = sigXX;                                                % stress storage
    sig(ilow:iup,2) = sigYY;
    sig(ilow:iup,3) = sigXY;
    
    strain(ilow:iup,1) = epsXX;                                             % strain storage
    strain(ilow:iup,2) = epsYY;
    strain(ilow:iup,3) = epsXY;
    
    ilow=ilow+nelblo;                                                       % final block adjustment
    if(ib==nbloc-1)
        nelblo=nels-iup;
        invJx=zeros(nelblo, nD);
        invJy=zeros(nelblo, nD);
    end
    iup=iup+nelblo;

end
uv = [sum(uvw(etpl*nD-1),2) sum(uvw(etpl*nD-0),2)]/nen;                     % Gauss point positions
xy = zeros(nels,nD);
for i = 1:nen
  xy = xy+[coord(etpl(:,i),1) coord(etpl(:,i),2)]/nen;                      % displaced Gauss point positions
end
fprintf('%s\n','   stresses, strains and displacements determined');