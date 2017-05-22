function [smpl, acceptance, energy, deformationV, maxTheta] = ...
    annealing(opt, u0, nsamples)
% Simulated annealing, with temperature descreasing as 1/iter
% Use MH_PLOT to generate plots.
% 
% last modified on Feb 24, 2017
% yun
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% opt - options
% u0 - starting point
% nsamples - total number of samples to take
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% smple - samples taken
% acceptance - acceptance rate at each step
% energy - energy of unit cell at each step
% deformationV - matrix containing vectors of deformation in each sample
% maxTheta - the maximum of deformation of hinge angles in each iteration
% 
% and also dumps all data in csv files
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


[unitCell,extrudedUnitCell,opt]=buildGeometry(opt);
[T]=linearConstr(unitCell,extrudedUnitCell,opt);
delta = opt.delta;
pdf = @(b, e) exp(-b * e);

% initialising...
smpl = zeros(nsamples, size(T,2));
acceptance = zeros(nsamples, 1);
energy = zeros(nsamples, 1);
deformationV = zeros(nsamples, size(T,1));
maxTheta = zeros(nsamples, 1);

u = u0(:);
% main loop
for ct = 1:nsamples
    [e,~,~,~,theta] = Energy(u, T, extrudedUnitCell, opt);
    beta = opt.beta * (1 - (ct-1)/nsamples); % gradually decrease temperature
    candidate = u + randn(size(T,2), 1) * delta;
    dTheta = theta - extrudedUnitCell.theta;
    % make sure that no faces are crossing over each other
    while ~isempty(find(abs(theta)>pi, 1))
        candidate = u + randn(size(T,2), 1) * delta;
    end
    
    [can_e,~,~,~,can_theta] = Energy(candidate, T, extrudedUnitCell, opt);
    acceptance(ct) = min(1, pdf(beta, can_e)/pdf(beta, e));
    
    if rand < acceptance(ct)
        u = candidate;
        e = can_e;
        dTheta = can_theta - extrudedUnitCell.theta;
    end
    
    smpl(ct,:) = u';
    energy(ct) = e;
    V = T * u;
    deformationV(ct,:) = V;
    [~,maxIdx] = max(abs(dTheta));
    maxTheta(ct) = dTheta(maxIdx);
end

if strcmp(opt.saveCSV, 'on')
    header = strcat('sim_anneal() ', ' Kface=', num2str(opt.Kface), ...
        ' Khinge=', num2str(opt.Khinge), ' Kedge=', num2str(opt.Kedge));
    saveCSV(opt, header, 'sim_anneal', 'smpl', smpl, ...
        'acceptance', acceptance, 'energy', energy, ...
        'deformationV', deformationV, 'maxTheta', maxTheta);
end

%% example annealing
% clear; close all
% opt=initOpt('inputType','individual',...
%             'template','truncated tetrahedron',...
%             'plot','result',... %'scale', 1,... % deformation
%             'saveMovie', 'off', 'safeMovieAntiAlias', 0, ...
%             'saveCSV', 'off', 'saveGraph', 'off', ... graph from mh_plot
%             'beta', .1, 'delta', .005, ...
%             'interval',1,'saveFig','off','periodic','off',... 
%             'Khinge',0,'Kedge',100,'Kface',100,'KtargetAngle',1,...
%             'date', datestr(now, 'mmm-dd-yyyy'),...
%             'time', datestr(now,'HH-MM-SS'));
% [unitCell,extrudedUnitCell,opt]=buildGeometry(opt);
% [T]=linearConstr(unitCell,extrudedUnitCell,opt);
% u0 = zeros(1, size(T,2)); % starting point
% nsamples = 1e2;
% [smpl, acceptance, energy, deformationV, maxTheta] = ...
%     annealing(opt, u0, nsamples);
% mh_plot(nsamples, energy, deformationV, opt, acceptance, maxTheta)