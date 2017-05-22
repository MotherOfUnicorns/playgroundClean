function [acceptance, energy, deformationV, maxTheta] = ...
    annealing(opt, u0, nsamples)
% [acceptance,energy,deformationV maxTheta] = annealing(opt,u0,nsamples)
% 
% Simulated annealing, with temperature descreasing as 1/iter
% Use MH_PLOT to generate plots.
% Also dumps all data in csv files
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% opt      - options
% u0       - starting point
% nsamples - total number of samples to take
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% acceptance   - acceptance rate at each step
% energy       - energy of unit cell at each step
% deformationV - matrix containing vectors of deformation in each sample
% maxTheta     - the maximum of deformation of hinge angles in each 
%                iteration
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% EXAMPLE
% clear; close all
% opt=initOpt('inputType','individual',...
%             'template','truncated tetrahedron',...
%             'plot','result',... %'scale', 1,... % deformation
%             'saveMovie', 'off', 'safeMovieAntiAlias', 0, ...
%             'saveCSV', 'off', 'saveGraph', 'off', ... graph from mh_plot
%             'beta', .1, 'delta', .005, ...
%             'interval',1,'saveFig','off','periodic','off',... 
%             'Khinge',0,'Kedge',100,'Kface',100,'KtargetAngle',1,...
%             'constrEdge','off', 'constrFace','off', ...
%             'date', datestr(now, 'mmm-dd-yyyy'),...
%             'time', datestr(now,'HH-MM-SS'));
% [~,extrudedUnitCell,opt]=buildGeometry(opt);
% u0 = zeros(1, 3 * size(extrudedUnitCell.node, 1)); % starting point
% nsamples = 1e2;
% [acceptance, energy, deformationV, maxTheta] = ...
%     annealing(opt, u0, nsamples);
% mh_plot(nsamples, energy, deformationV, opt, acceptance, maxTheta)
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on May 22, 2017
% yun


% probability density distribution is proportional to this
pdf = @(b, e) exp(-b * e);

% initialising...
delta = opt.delta;
[~,extrudedUnitCell,opt]=buildGeometry(opt);
N = 3 * size(extrudedUnitCell.node, 1); % number of variables
acceptance = zeros(nsamples, 1);
energy = zeros(nsamples, 1);
deformationV = zeros(nsamples, N);
maxTheta = zeros(nsamples, 1);

u = u0(:);
% main loop
for ct = 1:nsamples
    [e,~,~,~,~,theta] = Energy(u, extrudedUnitCell, opt);
    beta = opt.beta * (1-(ct-1)/nsamples); % gradually decrease temperature
    candidate = u + randn(N, 1) * delta;
    dTheta = theta - extrudedUnitCell.theta;
    % make sure that no faces are crossing over each other
    while ~isempty(find(abs(theta)>pi, 1))
        candidate = u + randn(N, 1) * delta;
    end
    
    [can_e,~,~,~,~,can_theta] = Energy(candidate,extrudedUnitCell,opt);
    acceptance(ct) = min(1, pdf(beta, can_e)/pdf(beta, e));
    
    if rand < acceptance(ct)
        u = candidate;
        e = can_e;
        dTheta = can_theta - extrudedUnitCell.theta;
    end
    

    energy(ct) = e;
    deformationV(ct,:) = u';
    [~,maxIdx] = max(abs(dTheta));
    maxTheta(ct) = dTheta(maxIdx);
end


% dump csv file
if strcmp(opt.saveCSV, 'on')
    header = strcat('sim_anneal() ', ' Kface=', num2str(opt.Kface), ...
        ' Khinge=', num2str(opt.Khinge), ' Kedge=', num2str(opt.Kedge));
    saveCSV(opt, header, 'sim_anneal', ...
        'acceptance', acceptance, 'energy', energy, ...
        'deformationV', deformationV, 'maxTheta', maxTheta);
end