strain = [ 0.0517    0.0394    0.0759    0.0660;
          -0.2603   -0.1958   -0.3666   -0.3428;
           0.0094    0.0037   -0.0157   -0.0086 ];
P_arg = -3;
varargin = '';

% Forward problem for 2x2 mesh:
%--------------------------------------------------------------------------
% Author: Stefan Szyniszewski and Edward Street
% Date:   21/06/2024
% Description: the script to solve the inverse problem, i.e. iteratively 
% computing Young's modulus (E) and Poisson's ratio (nu) for known
% boundary conditions, loads, and with given element strain values.
%
% *Requires the forward_problem_2x2_elems_ES.m function in the same path
%
% Varargin:
%         'notime'  = omits computing the run time
%
% Example strain values for E = [25 35 15 15], nu = [0.2 0.2 0.2 0.2], and P = -3:
%
%                     elem.1    elem.2    elem.3    elem.4
%                        |         |         |         |
%         strain =  [ 0.0517    0.0394    0.0759    0.0660;  % - ε_xx
%                    -0.2603   -0.1958   -0.3666   -0.3428;  % - ε_yy
%                     0.0094    0.0037   -0.0157   -0.0086 ] % - γ_xy
%
%--------------------------------------------------------------------------
% clearvars; close all; 
% clc;
% [thisPath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
% cd(thisPath);                                                             % change directory to current path
% addpath('functions');                                                       % add path to functions folder

%% Start clock
if ~any(strcmp(varargin,'notime'))
    tic;
end

%% Compute principal strains for each element in the measured solution
epsilon1 = max((strain(1,:)+strain(2,:))/2+sqrt(((strain(1,:)-strain(2,:))/2).^2+(strain(3,:)/2).^2), ...
               (strain(1,:)+strain(2,:))/2-sqrt(((strain(1,:)-strain(2,:))/2).^2+(strain(3,:)/2).^2));
epsilon2 = min((strain(1,:)+strain(2,:))/2+sqrt(((strain(1,:)-strain(2,:))/2).^2+(strain(3,:)/2).^2), ...
               (strain(1,:)+strain(2,:))/2-sqrt(((strain(1,:)-strain(2,:))/2).^2+(strain(3,:)/2).^2));

%% Initial guess for E, nu in all elements (i = 0)
nu0 = 0.4;
E0 = 10;
i = 0; % iteration counter

%% Get first solution using initial guess
% Compute forward solution
[~,epsilon,sigma] = forward_problem_2x2_elems_ES([E0 E0 E0 E0],[nu0 nu0 nu0 nu0],P_arg,'noprint','noplot','notime');

%% Compute the L-infinity norm between the measured and computed strains
change = max(abs(epsilon-strain),[],'all');
fprintf('i = %3d:\t E = [%.2f %.2f %.2f %.2f],\t nu = [%.2f %.2f %.2f %.2f],\t L_\x221E = %f\n',i,E0,E0,E0,E0,nu0,nu0,nu0,nu0,change)

%% Compute principal stresses in each element in the first solution
sigma1 = max((sigma(1,:)+sigma(2,:))/2+sqrt(((sigma(1,:)-sigma(2,:))/2).^2+(sigma(3,:)).^2), ...
             (sigma(1,:)+sigma(2,:))/2-sqrt(((sigma(1,:)-sigma(2,:))/2).^2+(sigma(3,:)).^2));
sigma2 = min((sigma(1,:)+sigma(2,:))/2+sqrt(((sigma(1,:)-sigma(2,:))/2).^2+(sigma(3,:)).^2), ...
             (sigma(1,:)+sigma(2,:))/2-sqrt(((sigma(1,:)-sigma(2,:))/2).^2+(sigma(3,:)).^2));

%% Update E and nu
% disp(epsilon1)
% disp(epsilon2)
% disp(sigma1)
% disp(sigma2)
nu_new = (sigma1.*epsilon2-sigma2.*epsilon1)./(sigma2.*epsilon2-sigma1.*epsilon1);
nu_new = max(-0.999,min(0.4999,nu_new));
E_new = (sigma1-nu_new.*epsilon2)./epsilon1;
E_new = max(0.001,E_new);

%% Repeat the process until the solution strain is close enough to the measured strain
% setting a tolerance for the change value
tol = 1e-3;
Linf = nan(1,100);
while change>tol && i<100
    i = i+1;
    % Compute forward solution
    [~,epsilon,sigma] = forward_problem_2x2_elems_ES(E_new,nu_new,P_arg,'noprint','noplot','notime');
    % Compute the L-infinity norm between the measured and computed strains
    change = max(abs(epsilon-strain),[],'all');
    Linf(i) = change; % record the L_infinity norm to plot convergence
    fprintf('i = %3d:\t E = [%.2f %.2f %.2f %.2f],\t nu = [%.2f %.2f %.2f %.2f],\t L_\x221E = %f\n',i,E_new,nu_new,change)
    % Compute principal stresses in each element in the first solution
    sigma1 = max((sigma(1,:)+sigma(2,:))/2+sqrt(((sigma(1,:)-sigma(2,:))/2).^2+(sigma(3,:)/2).^2), ...
                 (sigma(1,:)+sigma(2,:))/2-sqrt(((sigma(1,:)-sigma(2,:))/2).^2+(sigma(3,:)/2).^2));
    sigma2 = min((sigma(1,:)+sigma(2,:))/2+sqrt(((sigma(1,:)-sigma(2,:))/2).^2+(sigma(3,:)).^2), ...
                 (sigma(1,:)+sigma(2,:))/2-sqrt(((sigma(1,:)-sigma(2,:))/2).^2+(sigma(3,:)).^2));
    % Update E and nu 
    nu_new = (sigma1.*epsilon2-sigma2.*epsilon1)./(sigma2.*epsilon2-sigma1.*epsilon1);
    nu_new = max(-0.999,min(0.4999,nu_new));
    E_new = (sigma1-nu_new.*epsilon2)./epsilon1;
    E_new = max(0.001,E_new);
end
% Plot final solution deformation
[~,~,~] = forward_problem_2x2_elems_ES(E_new,nu_new,P_arg,'noprint','notime');

%% Return final E and nu values
E_return = E_new;
nu_return = nu_new;

%% Plot convergence
figure;
plot(1:i,Linf(1:i))
xlabel('i')
ylabel('L_infinity Norm')
title('Convergence')

%% Stop clock
if ~any(strcmp(varargin,'notime'))
    toc;
end