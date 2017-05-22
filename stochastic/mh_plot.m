function mh_plot(nsamples, energy, deformationV, opt, acceptance, maxTheta)
% mh_plot(nsamples, energy, deformationV, opt, acceptance, maxTheta)
% 
% Plots the results from MH.m, breating.m, annealing.m and many others
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% nsamples     - number of samples/iterations
% energy       - an array containing the energy at each step
% deformationV - a matrix containing all deformations
% opt          - options
% acceptance   - an array containing the acceptance rate at each step
%                (optional)
% maxTheta     - maximum deformation in angles of each hinge at each step
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% Produces an animation of the extruded unit cell, as well as a graph of
% all the input values.
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% EXAMPLE
% opt=initOpt('inputType','individual',...
%             'template','truncated tetrahedron',...
%             'plot','result',... %'scale', 1,... % deformation
%             'saveMovie', 'on', 'safeMovieAntiAlias', 0, ...
%             'saveCSV', 'on', 'saveGraph', 'on', ... graph from mh_plot
%             'beta', 1, 'delta', .02, ...
%             'interval',1,'saveFig','off','periodic','off',... 
%             'Khinge',0,'Kedge',1,'Kface',100,'KtargetAngle',1,...
%             'constrEdge','off', 'constrFace','off', ...
%             'date', datestr(now, 'mmm-dd-yyyy'),...
%             'time', datestr(now,'HH-MM-SS'));
% [~,extrudedUnitCell,opt]=buildGeometry(opt);
% u0 = zeros(1, 3*size(extrudedUnitCell.node, 1)); % starting point
% nsamples = uint32(2e2);
% tic;
% [acceptance, energy, deformationV, maxTheta] = mh(opt, u0, nsamples);
% toc;
% max(maxTheta)
% tic;
% mh_plot(nsamples, energy, deformationV, opt, acceptance, maxTheta)
% toc;
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on May 22, 2017
% yun

% input validation...
if ~exist('acceptance', 'var')
    has_acceptance = false;
    acceptance = [];
else
    has_acceptance = true;
end
if ~exist('maxTheta', 'var')
    has_maxTheta = false;
    maxTheta = [];
else
    has_maxTheta = true;
end

% setting up...
[unitCell,extrudedUnitCell,opt]=buildGeometry(opt);
date = opt.date;
time = opt.time;


deformation = zeros(nsamples, 1);
for ct = 1:length(deformation)
    deformation(ct) = norm(deformationV(ct,:));
end


h = figure(1); hold on;
%title(strcat('beta = ', num2str(beta), ...
%    ', nsamples = ', num2str(nsamples)));

if has_acceptance
    subplot(3,2,1)
    plot(acceptance, '.');
    xlabel('n')
    ylabel('acceptance')
end


if has_maxTheta
    subplot(3,2,2)
    plot(maxTheta, '.');
    xlabel('n')
    ylabel('maxTheta')
    
    % also indicate the largest 10 deformation
    [~,sortIndex] = sort(abs(maxTheta),'descend');
    maxIndex = sortIndex(1:5);
    hold on;
    plot(maxIndex, maxTheta(maxIndex), 'xr')
    hold off;
end

subplot(3,2,3)
plot(energy, '.');
xlabel('n')
ylabel('energy')

subplot(3,2,4)
histogram(energy, 'orientation', 'horizontal',...
    'NumBins', max(5, ceil(nsamples/100)));
xlabel('count')
ylabel('energy')

subplot(3,2,5)
plot(deformation,'.');
xlabel('n')
ylabel('deformation')

subplot(3,2,6)
plot(deformation, energy, '.');
ylabel('energy')
xlabel('deformation')

% creates folder if doesn't exist
nameFolder = [pwd, '/Results/', date];
if exist(nameFolder, 'dir')==0
    mkdir(nameFolder);
end

% save the graph...
if strcmp(opt.saveGraph, 'on')
    savefig(h, [nameFolder '/', ...
        strcat(opt.template, '_', time, '_', 'beta=', num2str(opt.beta),...
        'delta=', num2str(opt.delta),'.fig')]);
end

% only update animation every 10 iterations
updateRate = 10;
timesteps = 1:updateRate:nsamples;

result.numMode = length(timesteps);
for ct = 1:result.numMode
    idx = (ct-1) * updateRate + 1;
    result.deform(ct).V = [deformationV(idx, 1:3:end)' ...
        deformationV(idx, 2:3:end)' deformationV(idx, 3:3:end)'];
    result.deform(ct).Ve = deformationV(idx,:)';
end

% also include the last frame when nsamples is not a multiple of updateRate
if mod(nsamples-1, updateRate) ~= 0
    timesteps = [timesteps, nsamples];
    result.numMode = result.numMode + 1;
    result.deform(result.numMode).V = [deformationV(end, 1:3:end)' ...
        deformationV(end, 2:3:end)' deformationV(end, 3:3:end)'];
    result.deform(result.numMode).Ve = deformationV(end,:)';
end

outputResults(unitCell,extrudedUnitCell,result,timesteps,opt);