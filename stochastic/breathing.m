function [acceptance, energy, deformationV, maxTheta] = ...
    breathing(opt, u0, nsamples, nrepeats)
% [acceptance, energy, deformationV, maxTheta] = ...
%    breathing(opt, u0, nsamples, nrepeats)
%
% Metropolis-Hastings with alternating high and low temperatures
% Use MH_PLOT to generate plots.
% Also dumps all data in csv files.
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% opt - options
% u0 - starting point
% nsamples - total number of samples to take
% nrepeats - how many cycles to repeat (one cycle consists of one period of
%       high temperature and one period of low temperature)
% ***** the total number of iterations will be (nsamples * nrepeats) *****
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% smple - samples taken
% acceptance - acceptance rate at each step
% energy - energy of unit cell at each step
% deformationV - matrix containing vectors of deformation in each sample
% maxTheta - the maximum of deformation of hinge angles in each iteration
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% EXAMPLE
% clear; close all
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
% u0 = zeros(1, 3*size(extrudedUnitCell.node,1)); % starting point
% nsamples = uint32(5e1);
% nrepeats = 3;
% tic;
% [acceptance, energy, deformationV, maxTheta] = ...
%     breathing(opt, u0, nsamples, nrepeats);
% toc;
% mh_plot(nsamples*nrepeats, energy, deformationV, opt, acceptance, maxTheta)
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on May 22, 2017
% yun

% input validation and initialising...
if ~exist('nrepeats', 'var')
    nrepeats = 1;
end
temp_factor = 5; % how much beta is changed in the low-temperature round
tot_iter = nsamples * nrepeats;

[~,extrudedUnitCell,opt]=buildGeometry(opt);
N = 3 * size(extrudedUnitCell.node, 1); % number of variables
acceptance = zeros(tot_iter, 1);
energy = zeros(tot_iter, 1);
deformationV = zeros(tot_iter, N);
maxTheta = zeros(tot_iter, 1);

for cycle = 1:nrepeats
    % first round with a higher temperature
    nsamples1 = ceil(nsamples / 2);
    [acceptance1,energy1,deformationV1,maxTheta1] = mh(opt,u0,nsamples1);

    % second round with a lower temperature
    u0 = deformationV1(end,:);
    opt.beta = opt.beta * temp_factor;
    nsamples2 = nsamples - nsamples1;
    [acceptance2,energy2,deformationV2,maxTheta2] = mh(opt,u0,nsamples2);

    % restoring beta for the (possible) next iteration
    opt.beta = opt.beta / temp_factor;

    % saving array
    idx0 = (cycle-1) * nsamples + 1; 
    idx1 = cycle * nsamples;
    acceptance(idx0:idx1)=[acceptance1;acceptance2];
    energy(idx0:idx1)=[energy1;energy2];
    deformationV(idx0:idx1,:)=[deformationV1;deformationV2];
    maxTheta(idx0:idx1) = [maxTheta1; maxTheta2];
end
    

if strcmp(opt.saveCSV, 'on')
    header = strcat('annealing() ', ' Kface=', num2str(opt.Kface), ...
        ' Khinge=', num2str(opt.Khinge), ' Kedge=', num2str(opt.Kedge));
    saveCSV(opt, header, 'annealing', 'acceptance', acceptance,...
        'energy', energy, 'deformationV', deformationV);
end