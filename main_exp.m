%plane stress field and elastic property solver
%--------------------------------------------------------------------------
% Author: Stefan Szyniszewski
% Date:   06/07/2021
% Description: Script to computer stress field and elastic constants
% in a material sample, based on a known (measured) strain field, epsA.
% The script iterates until it converges within a given tolerance or
% reaches a maximum number of iteraitons (default set at 50).
%
%--------------------------------------------------------------------------
% MAIN
%--------------------------------------------------------------------------
% Input(s):
%   - Strain.xlsx   - file with the known strain information
%--------------------------------------------------------------------------
% Ouput(s):
%   - .\results\Results.xlsx - Excell file with numeric results
%   - .\results\epsilon-xx - input.png  - images with results
%--------------------------------------------------------------------------
% See also:
% setupmesh         - problem set up
% LEfe              - finite element solver
% updateElasticProp - elastic property update
%--------------------------------------------------------------------------
clearvars; close all;
clc;
% [thisPath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
% cd(thisPath);                                                             % change directory to current path
addpath('functions');                                                       % add path to functions folder
clear; tic;                                                                 % clear all breakpoints and start stopwatch for computatotional efficiency calcs

%% Physical problem description:
% *****************************
pressure = -0.06;                                                           % pressure in [Pa]
element_size = 0.8592;                                                      % physical length (of pixel) in [m] or other consistent units
bc_type = 'pressure';                                                       % type of boundaries at the bottom: 'pressure' or 'roller' or 'fixed'
   
%% Read strain inputs:
%****************************
% readCSVfile
epsA = xlsread('data\exp_trial\trial_15.csv');                              % strain input: xx, yy, xy
epsA = epsA(:,3:end);
% reshape_the_data;
% Check the data

%% Set solver tolerances:
% **********************
itMax = 18;                                                                 % maximum number of iterations
tol   = 1e-7;                                                               % tolerance
filter_type = 'None';                                                       % 'None','Gaussian' or 'MovingAverage'
filter_size = 3;                                                            % 3, 5, or any higher odd number.

%% Set up the physical model (based on the user input and characteristics of the supplied data):
% *************************************************************************
nels_x = 47;                                                                % number of elements in the x direction
nels_y = 91;                                                                % number of elements in the y direction
% **************
nels = nels_x * nels_y;                                                     % total number of elements
l_x = nels_x * element_size;                                                % domain size in x-direction
l_y = nels_y * element_size;                                                % domain size in y-direction

%      Y
%      ^
%      |     pressure load
%      |     -------------
%      | || || || || || || || ||
%      | \/ \/ \/ \/ \/ \/ \/ \/        |
%      |------------------------.       -
%      | 26 | 27 | 28 | 29 | 30 |       |
%      |------------------------|       |
%      | 21 | 22 | 23 | 24 | 25 |       |
%      |------------------------|       |
%      | 16 | 17 | 18 | 19 | 20 |      l_x   (nels_y)
%      |------------------------|       |
%      | 11 | 12 | 13 | 14 | 15 |       |
%      |------------------------|       |
%      |  6 |  7 |  8 |  9 | 10 |       |
%      |------------------------|       |
%      |  1 |  2 |  3 |  4 |  5 |       |
%      -----------------------------------------> X
%      ^    ^    ^    ^    ^    ^
%       simply supported boundary
%
%     -/-------- l_y -----------/--
%    (nels_x - number of elements in X-direction)

% Compute average Young modulus and Poisson ratio as starting points:
eps_xx_avg = mean(epsA(:,1));                                               % average tranverse strain
eps_yy_avg = mean(epsA(:,2));                                               % average axial (loading direction) strain

Lx = nels_x * element_size;                                                 % length of the domain
E_init = -pressure / eps_yy_avg;                                            % average Young modulus
v_init = -eps_xx_avg/eps_yy_avg;                                            % average Poisson ratio

%% Set FE parameters based on the user input and characteristics of the supplied data:
mesh = setupmesh(l_x,l_y,nels_x,nels_y,pressure,E_init,v_init,bc_type);     % initialize mesh

%% Set data settings:
% *******************
% switch for saving incremental, iterative data and plots (for debugging
% only),

% Add strain convergence norm plots to the data for eps_xx, and eps_yy.

% create theia_log.txt file with informations about the progress of
% iterations, e.g. iter_num = ..., eps_xx_resid_norm = ...,
% eps_yy_resid_norm = ...

%% Material iterations
error = 2*tol;                                                              % initial error
delta_error = -tol;                                                         % initialize delta_error as negative number to start while loop
itnum = 0;                                                                  % zero iteration counter

while error > tol && itnum < itMax && delta_error < 0                       % while loop (as long as error is decreasing)
    itnum = itnum+1;                                                        % iteration counter
    fprintf('\n%s%8i\n','   iteration number     ',itnum);                  % print iteration number
    [sig,epsH,~] = LEfe(mesh);                                              % finite element solver
    for i=1:3                                                               % Gaussian smoothing of experimental noise
        if (strcmp(filter_type,'None') ~= 1)
            sig(:,i) = smoothVector(sig(:,i),nels_y,nels_x,filter_type,filter_size);
            epsH(:,i) = smoothVector(epsH(:,i),nels_y,nels_x,filter_type,filter_size);
        end
    end
    [E,v] = updateElasticProp(sig,epsA);                                    % solve for elastic properties
    mesh.E = E;                                                             % update Young's modulus
    mesh.v = v;                                                             % update Poisson's ratio
    error_old = error;                                                      % error from previous iteration
    error = norm(epsH(:,1:2)-epsA(:,1:2))/norm(epsA(:,1:2));                % normalised strain error
    if itnum == 1
        delta_error = -2*tol;                                               % keep delta_error negative to keep the while loop going
    else
        delta_error = error - error_old;                                    % check if delta_error < 0, which means the error is decreasing
    end
    fprintf('%s%8.3e\n','   error                   ',error);               % return strain error
    figure(1001);                                                           % figure number
    matrix = vector2matrix(sig(:,2),nels_y,nels_x);                         % convert vector to matrix for plotting
    imagesc(matrix);colorbar; colormap; axis equal; axis off;               % plot color map
    drawnow limitrate nocallbacks;                                          % drawnow with 20 frames per second limit
    title( strcat('\sigma - yy - iteration:',num2str(itnum)) );
    errorit(itnum)=error;
end

%% Save final results to Excel file
mydir  = pwd;                                                               % current working directory
output_dir = strcat(mydir,'\results\');                                     % create iterations folder
if (isfolder(output_dir))                                                   % do nothing, the folder exists
else
    mkdir(output_dir);                                                      % create direcctory if it does not exist
end
% eps_xx,eps_yy,eps_yy,sig_xx,sig_yy,sig_xy,E,v,eps_xx_H,eps_yy_H,eps_xy_H
% ------- input ------, ----- stress ------,mat,---- numerical strains ---
results = [epsA,sig,E,v,epsH];
xlswrite(strcat(output_dir,'results.xlsx'),results);                        % write results to Excell file

%% Plot the results
fig_num = 0;                                                                % initialize figure number

% Matched strains:
plot_titles = {'epsilon-xx - matched', 'epsilon-yy - matched', 'epsilon-xy - matched'};
fig_num = fig_num +1; %close(fig_num);                                      % figure number
my_subplot(fig_num,plot_titles,epsA,nels_x, nels_y);
saveas(gcf,strcat(output_dir,'strains-matched','.png'));                      % save images as png

% Strains - input:
plot_titles = {'epsilon-xx - input', 'epsilon-yy - input', 'epsilon-xy - input'};
fig_num = fig_num +1; %close(fig_num);                                      % figure number
my_subplot(fig_num,plot_titles,epsA,nels_x, nels_y);
saveas(gcf,strcat(output_dir,'strains-input','.png'));                      % save images as png

% Elastic constants:
fig_num = fig_num +1; %close(fig_num);                                      % figure number
sub = my_subplot(fig_num,{'Young', 'Poisson'},[E,v],nels_x, nels_y);
saveas(gcf,strcat(output_dir,'elastic_constants','.png'));                  % save images as png

% Stress:
fig_num = fig_num +1; %close(fig_num);                                      % figure number
my_subplot(fig_num,{'sigma-xx', 'sigma-yy', 'sigma-xy'},sig,nels_x, nels_y);
saveas(gcf,strcat(output_dir,'stress','.png'));                             % save images as png file

figure(fig_num+1)
hold on;
%scatter(1:1:length(errorit),errorit,200,'b','s','filled')
plot(1:1:length(errorit),errorit,'-s','MarkerSize',20,...
    'MarkerEdgeColor','b', 'MarkerFaceColor','b','linewidth',2,'color','k')
box on;
xlabel('Iteration')
ylabel('RVSE/VE_{avg}')
set(gca,'fontsize',32,'linewidth',2);
set(gcf,'position',[10,20,720,600])

%% FUNCTIONS
% plotting function:

function [sub] = my_subplot(fig_num,plot_titles,data,nels_x, nels_y)
    
    num_subplots = size(plot_titles,2);                                     % each gap is 20% of each plot
    size_subplot_x = (0.2 + num_subplots + (num_subplots-1)*0.2 + 0.2) * nels_x;
    size_subplot_y = (0.2 + 1.0 + 0.2) * nels_y;
    amplify_fig = num_subplots * 1/2;                                       % amplify the figure for larger number of subplots

    f = figure(fig_num);  close(fig_num); f = figure(fig_num);
    width_fig = f.Position(3); height_fig = f.Position(4);
    scale_fig_size_x = width_fig/size_subplot_x *amplify_fig;
    scale_fig_size_y = height_fig/size_subplot_y *amplify_fig;
    scale_fig_size = min(scale_fig_size_x,scale_fig_size_y);
%     f.Position(1) = 1/amplify_fig * f.Position(1);
    f.Position(2) = 1/num_subplots * f.Position(2);
    f.Position(3) = size_subplot_x * scale_fig_size_x;
    f.Position(4) = size_subplot_y * scale_fig_size_y;
    
    pos_subplot = zeros(num_subplots,4);                                    % placeholder for the position
    pos_subplot_1 = 0.06;                                                    % positions of the figures
    for ii=1:num_subplots
        pos_subplot(ii,:) = [ pos_subplot_1 0.15 nels_x/size_subplot_x nels_y/size_subplot_y];
        pos_subplot_1 = pos_subplot_1 + (1.2*nels_x)/size_subplot_x;
    end   

    for i=1:num_subplots                                                    % plot each subplot figure
        sub(i) = subplot(1,size(plot_titles,2),i); subplot('Position',pos_subplot(i,:))
        hold on;
        matrix = vector2matrix(data(:,i),nels_y,nels_x);                    % convert vector to matrix for plotting
        imagesc(matrix);colorbar; colormap(jet); %colormap(hsv);            % plot color map
        axis equal; axis off;
    %     caxis([minSigma maxSigma]);
        title( plot_titles{i} );                                            % add plot title
    %   title( strcat('\',plot_titles{i}) );                                % add plot title
    end
end
