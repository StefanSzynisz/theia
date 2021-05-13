function [N] = SF(gp,ngp,nen)

%Finite element Shape Functions (SF)
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   07/11/2019
% Description:
% Function to determine the shape functions of a four noded quarilaterial
% element at specific Gauss point locations. 
%
%--------------------------------------------------------------------------
% [N] = SF(gp,ngp,nen)
%--------------------------------------------------------------------------
% Input(s):
% gp     - Gauss point number of range required
% ngp    - number of Gauss points
% nen    - number of element nodes
%--------------------------------------------------------------------------
% Ouput(s);
% N      - shape function matrix
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------

tol=1e-3; 
if ngp==4
  g2=1/sqrt(3);
  xsi=[-1 -1 1 1].'*g2;
  eta=[-1 1 1 -1].'*g2;  
else
  xsi=0;
  eta=0;
end
N=zeros(nen,length(gp));
if abs(nen-4)<tol
  N(1,:)=(1-xsi(gp)).*(1-eta(gp))/4;
  N(2,:)=(1-xsi(gp)).*(1+eta(gp))/4;
  N(3,:)=(1+xsi(gp)).*(1+eta(gp))/4;
  N(4,:)=(1+xsi(gp)).*(1-eta(gp))/4;
end
