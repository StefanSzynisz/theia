% Estimation of mechanial properties of bones
%--------------------------------------------------------------------------
% Author: Stefan Szyniszewski
% Date:   23/12/2021
% Description: Estimation of mechanial properties of bones
%
%--------------------------------------------------------------------------
% MAIN
%--------------------------------------------------------------------------
% Input(s):
%   - 
%--------------------------------------------------------------------------
% Ouput(s):
%   - 
%--------------------------------------------------------------------------

E_bone = 9.0 * 10^9;                                                        % Pa, Modul Younga kosci
F = 6.0 * 9.81;                                                             % N, Applied force
A = 3.14 * (0.02)^2/4;                                                      % m^2, Approximate area of the bone cross-section

sig = F/A                                                                   % Applied stress
eps = sig/E_bone                                                            % Corresponding strain


