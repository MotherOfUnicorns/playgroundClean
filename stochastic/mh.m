function [acceptance, energy, deformationV, maxTheta] = ...
    mh(opt, u0, nsamples)
% [acceptance, energy, deformationV, maxTheta] = mh(opt, u0, nsamples)
%
% An implementation of the Metropolis-Hastings algorithm, where the
% proposal distribution is fixed to be a normal distribution.

% Uses MH_PLOT to generate plots,
% and also dumps all data in csv files.
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
% maxTheta     - maximum deformation of hinge angles in each iteration
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

beta = opt.beta;
delta = opt.delta;
[~,extrudedUnitCell,opt]=buildGeometry(opt);
pdf = @(e) exp(-beta * e);

%% initialising...
N = 3 * size(extrudedUnitCell.node, 1);
energy = zeros(nsamples, 1);
acceptance = zeros(nsamples, 1);
deformation = zeros(nsamples, 1);
deformationV = zeros(nsamples, N);
maxTheta = zeros(nsamples, 1);

%% main loop
u = u0(:);
for ct = 1:nsamples
    [e,~,~,~,~,theta] = Energy(u, extrudedUnitCell, opt);
    candidate = u + randn(N, 1) * delta;
    dTheta = theta - extrudedUnitCell.theta;
    % make sure that no faces are crossing over each other
    while ~isempty(find(abs(theta)>pi, 1))
        candidate = u + randn(N, 1) * delta;
    end
    
    [can_e,~,~,~,~,can_theta] = Energy(candidate,extrudedUnitCell,opt);
    acceptance(ct) = min(1, pdf(can_e)/pdf(e));
    
    % accept new candidate
    if rand < acceptance(ct)
        u = candidate;
        e = can_e;
        dTheta = can_theta - extrudedUnitCell.theta;
    end

    energy(ct) = e;

    deformationV(ct,:) = u';
    deformation(ct) = norm(u);
    [~,maxIdx] = max(abs(dTheta));
    maxTheta(ct) = dTheta(maxIdx);
end

if strcmp(opt.saveCSV, 'on')
    header = strcat('mh() beta=', num2str(opt.beta), ...
        ' delta=', num2str(opt.delta), ' Kface=', num2str(opt.Kface), ...
        ' Khinge=', num2str(opt.Khinge), ' Kedge=', num2str(opt.Kedge));
    saveCSV(opt, header, 'mh', 'acceptance', acceptance,...
        'energy',energy, 'deformationV',deformationV, 'maxTheta',maxTheta);
end