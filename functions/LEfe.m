function [sig,strain,uv] = LEfe(mesh,itnum)

%LEfe - Linear Elastic finite element solver
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   02/04/2020
% Description:
% Finite element solver for 2D plane stress analysis based on a vectorised
% implementation. 
%
%--------------------------------------------------------------------------
% [sig,strain] = LEFE(mesh,itnum)
%--------------------------------------------------------------------------
% Input(s):
% mesh   - structured array of mesh data
% itnum  - iteration number 
%--------------------------------------------------------------------------
% Ouput(s);
% sig     - Cauchy stress at the centre of each element (nels,3)
% strain  - strain at the centre of each element (nels,3)
% 
%--------------------------------------------------------------------------
% See also:
% 
% FORMKVEC            - global stiffness formation (vectorised)
% AVBCSVEC            - average boundary conditions (vectorised)
% LINEARSOLVER        - linear solver
% POSTPRO             - Cauchy stress determination
% MAKEVTK             - VTK mesh file generation
% MAKEVTKMP           - VTK stress data file generation
%--------------------------------------------------------------------------

iterationPlotInterval = 10;                                                 % plot every e.g. 5 iterations, 5, 10, 15, ..
numberOfDetailedIteration = 3;                                              % plot every iteration until this limit, eg. 1, 2,..5 

K    = formKvec(mesh);                                                      % vectorised stiffness calculation
fext = mesh.fext;                                                           % external force vector

uvw = linearSolver(K,mesh.bc,fext);                                         % displacement solution
[sig,strain,uv,xy] = postPro(mesh,uvw);                                     % vectorised stress/displacement calculation

mydir  = pwd;                                                               % current working directory
outputPath = strcat(mydir,'\iterations\');                                  % create iterations folder

coordinate_x = xy(:,1); coordinate_y = xy(:,2);

sigxx = sig(:,1);
min_sigxx = min(sigxx); max_sigxx = max(sigxx);
resolution = 0.01;
plotTitle = sprintf('Stress in the X direction-iteration %i',itnum);
figure(1); left_pos = 3; bott_pos = 12; %cm
% Add axis option (manual tuning):
saveDataPath = strcat(outputPath,'\sigxx\');
if itnum==1
    mkdir(saveDataPath);  
end

% Plot only 1,2, ..5, 15, 20, 25, ..
if (itnum <= numberOfDetailedIteration) || (rem(itnum,iterationPlotInterval)==0)
    surfPlot(coordinate_x,coordinate_y,sigxx,plotTitle,'', saveDataPath,resolution,...
        left_pos,bott_pos, min_sigxx, max_sigxx, 'colorbar_on');
end

sigyy = sig(:,2);
min_sigyy = min(sigyy); max_sigyy = max(sigyy);
resolution = 0.01;
plotTitle = sprintf('Stress in the Y direction-iteration %i',itnum);
figure(2); left_pos = 13; bott_pos = 12; %cm
% Add axis option (manual tuning):
saveDataPath = strcat(outputPath,'\sigyy\');
if itnum==1
    mkdir(saveDataPath);  
end

if (itnum <= numberOfDetailedIteration) || (rem(itnum,iterationPlotInterval)==0)
    surfPlot(coordinate_x,coordinate_y,sigyy,plotTitle,'', saveDataPath,resolution,...
        left_pos,bott_pos, min_sigyy, max_sigyy, 'colorbar_on');
end

sigxy = sig(:,3);
min_sigxy = min(sigxy); max_sigxy = max(sigxy);
resolution = 0.01;
plotTitle = sprintf('Shear stress-iteration %i',itnum);
figure(3); left_pos = 23; bott_pos = 12; %cm
% Add axis option (manual tuning):
saveDataPath = strcat(outputPath,'\sigxy\');
if itnum==1
    mkdir(saveDataPath);  
end

if (itnum <= numberOfDetailedIteration) || (rem(itnum,iterationPlotInterval)==0)
    surfPlot(coordinate_x,coordinate_y,sigxy,plotTitle,'', saveDataPath,resolution,...
        left_pos,bott_pos, min_sigxy, max_sigxy, 'colorbar_on');
end

strainxx = strain(:,1);
min_strainxx = min(strainxx); max_strainxx = max(strainxx);
resolution = 0.01;
plotTitle = sprintf('Strain in the X direction-iteration %i',itnum);
figure(4); left_pos = 3; bott_pos = 2; %cm
% Add axis option (manual tuning):
saveDataPath = strcat(outputPath,'\strainxx\');
if itnum==1
    mkdir(saveDataPath);  
end

if (itnum <= numberOfDetailedIteration) || (rem(itnum,iterationPlotInterval)==0)
    surfPlot(coordinate_x,coordinate_y,strainxx,plotTitle,'', saveDataPath,resolution,...
        left_pos,bott_pos, min_strainxx, max_strainxx, 'colorbar_on');
end

strainyy = strain(:,2);
min_strainyy = min(strainyy); max_strainyy = max(strainyy);
resolution = 0.01;
plotTitle = sprintf('Strain in the Y direction-iteration %i',itnum);
figure(5); left_pos = 13; bott_pos = 2; %cm
% Add axis option (manual tuning):
saveDataPath = strcat(outputPath,'\strainyy\');
if itnum==1
    mkdir(saveDataPath);  
end

if (itnum <= numberOfDetailedIteration) || (rem(itnum,iterationPlotInterval)==0)
    surfPlot(coordinate_x,coordinate_y,strainyy,plotTitle,'', saveDataPath,resolution,...
        left_pos,bott_pos, min_strainyy, max_strainyy, 'colorbar_on');
end

strainxy = strain(:,3);
min_strainxy = min(strainxy); max_strainxy = max(strainxy);
resolution = 0.01;
plotTitle = sprintf('Shear strain-iteration %i',itnum);
figure(6); left_pos = 23; bott_pos = 2; %cm
% Add axis option (manual tuning):
saveDataPath = strcat(outputPath,'\strainxy\');
if itnum==1
    mkdir(saveDataPath);  
end

if (itnum <= numberOfDetailedIteration) || (rem(itnum,iterationPlotInterval)==0)
    surfPlot(coordinate_x,coordinate_y,strainxy,plotTitle,'', saveDataPath,resolution,...
        left_pos,bott_pos, min_strainxy, max_strainxy, 'colorbar_on');
end

ux = uv(:,1);
min_ux = min(ux); max_ux = max(ux);
resolution = 0.01;
plotTitle = sprintf('Displacement in the X direction-iteration %i',itnum);
figure(7); left_pos = 3; bott_pos = -8; %cm
% Add axis option (manual tuning):
saveDataPath = strcat(outputPath,'\ux\');
if itnum==1
    mkdir(saveDataPath);  
end

if (itnum <= numberOfDetailedIteration) || (rem(itnum,iterationPlotInterval)==0)
    surfPlot(coordinate_x,coordinate_y,ux,plotTitle,'', saveDataPath,resolution,...
        left_pos,bott_pos, min_ux, max_ux, 'colorbar_on');
end

vy = uv(:,2);
min_vy = min(vy); max_vy = max(vy);
resolution = 0.01;
plotTitle = sprintf('Displacement in the Y direction-iteration %i',itnum);
figure(8); left_pos = 23; bott_pos = -8; %cm
% Add axis option (manual tuning):
saveDataPath = strcat(outputPath,'\vy\');
if itnum==1
    mkdir(saveDataPath);  
end

if (itnum <= numberOfDetailedIteration) || (rem(itnum,iterationPlotInterval)==0)
    surfPlot(coordinate_x,coordinate_y,vy,plotTitle,'', saveDataPath,resolution,...
        left_pos,bott_pos, min_vy, max_vy, 'colorbar_on');
end

E = mesh.E;
min_E = min(E); max_E = max(E);
resolution = 0.01;
plotTitle = sprintf('Youngs modulus-iteration %i',itnum);
figure(9); left_pos = 9; bott_pos = 12; %cm
% Add axis option (manual tuning):
saveDataPath = strcat(outputPath,'\E\');
if itnum==1
    mkdir(saveDataPath);  
end

if (itnum <= numberOfDetailedIteration) || (rem(itnum,iterationPlotInterval)==0)
    surfPlot(coordinate_x,coordinate_y,E,plotTitle,'', saveDataPath,resolution,...
        left_pos,bott_pos, min_E, max_E, 'colorbar_on');
end   

v = mesh.v;
min_v = min(v); max_v = max(v);
resolution = 0.01;
plotTitle = sprintf('Poisson ratio-iteration %i',itnum);
figure(10); left_pos = 19; bott_pos = 12; %cm
% Add axis option (manual tuning):
saveDataPath = strcat(outputPath,'\v\');
if itnum==1
    mkdir(saveDataPath);  
end

if (itnum <= numberOfDetailedIteration) || (rem(itnum,iterationPlotInterval)==0)
    surfPlot(coordinate_x,coordinate_y,v,plotTitle,'', saveDataPath,resolution,...
        left_pos,bott_pos, min_v, max_v, 'colorbar_on');
end

end