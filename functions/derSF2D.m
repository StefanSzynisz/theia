function [wp,dNr]=derSF2D(ngp,nen)

%Shapefunction derivatives and Gauss weights
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   06/05/2015
% Description:
% Function to provide the derivatives of the finite element shape functions
% with respect to the local coordinates within the element at the Gauss
% point positions and the weights associated with the Gauss points for 2D
% analyses.
%
%--------------------------------------------------------------------------
% [wp,dNr] = DERSF2D(ngp,nen)
%--------------------------------------------------------------------------
% Input(s):
% ngp   - number of Gauss points (total)
% nen   - number of element nodes
%--------------------------------------------------------------------------
% Ouput(s);
% wp    - Gauss point weights
% dNr   - local derivatives of the shape functions
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------

if ngp==1
    gp = [0 0];
    wp = 4;
elseif ngp==4
    gp = [-1 -1; -1 1; 1 1; 1 -1]/sqrt(3);
    wp = ones(4,1);
end
xsi=gp(:,1); eta=gp(:,2); r2=ngp*2;
if nen==4
    dNr(1:2:r2  ,1) = -1/4*(1-eta);
    dNr(1:2:r2  ,2) = -1/4*(1+eta);
    dNr(1:2:r2  ,3) =  1/4*(1+eta);
    dNr(1:2:r2  ,4) =  1/4*(1-eta);
    dNr(2:2:r2+1,1) = -1/4*(1-xsi);
    dNr(2:2:r2+1,2) =  1/4*(1-xsi);
    dNr(2:2:r2+1,3) =  1/4*(1+xsi);
    dNr(2:2:r2+1,4) = -1/4*(1+xsi);
end