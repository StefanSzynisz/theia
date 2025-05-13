function [E,v] = updateElasticProp(sig,epsA)

%Plane stress elastic properties update 
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   15/04/2020
% Description:
% Function to update the elastic material constants for a linear elastic
% material under a plane stress condition.  The function updates the
% Young's modulus and Poisson's ratio of the material based on the
% predicted stress and a reference strain state. 
%
%--------------------------------------------------------------------------
% [E,v] = UPDATEELASTICPROP(sig,epsA)
%--------------------------------------------------------------------------
% Input(s):
% sig   - predicted stress values (nels,3)
%   xx, yy and xy order ??
% epsA  - reference strain values (nels,3)
%   xx, yy and xy order ??
%--------------------------------------------------------------------------
% Ouput(s);
% E     - updated Young's modulus (nels,1)
% v     - updated Poisson's ratio (nels,1)
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------

nels = length(sig);                                                         % number of elements 
E    = zeros(nels,1);                                                       % zero Young's modulus
v    = zeros(nels,1);                                                       % zero Poisson's ratio

for nel = 1:nels                                                            % loop over elements
    sigM = [sig(nel,1) sig(nel,3); sig(nel,3) sig(nel,2)];                  % stress in matrix format
    s    = sort(eig(sigM),'descend');                                       % principal stress values
    epsM = [epsA(nel,1) epsA(nel,3)/2; epsA(nel,3)/2 epsA(nel,2)];          % strain in matrix format
    e    = sort(eig(epsM),'descend');                                       % principal strain values
    Poisson = (e(2)*s(1)-e(1)*s(2))/(e(2)*s(2)-e(1)*s(1));                  % updated Poisson's ratio
    Poisson = min(Poisson, 0.4999);
    Poisson = max(-0.999, Poisson);
    v(nel) = Poisson;                                                       % updated Poisson's ratio
    Young = (s(1)-v(nel)*s(2))/e(1);                                        % updated Young's modulus
    E(nel) = max(10^-3,Young);                                              % make sure Young modulus is greater than zero
end
fprintf('%s\n','   elastic properties updated');