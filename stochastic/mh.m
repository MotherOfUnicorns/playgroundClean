function [smpl, acceptance, energy, deformationV, maxTheta] = ...
    mh(opt, u0, nsamples)
% finding stable states of extruded unit cell with Metropolis-Hastings
% algorithm.
% Use MH_PLOT to generate plots.
% 
% last modified on Feb 24, 2017
% yun
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% opt - options
% u0 - starting point
% nsamples - total number of samples to take
% TO ADD: (prop) - type of proposal function
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


beta = opt.beta;
delta = opt.delta;
[unitCell,extrudedUnitCell,opt]=buildGeometry(opt);
[T]=linearConstr(unitCell,extrudedUnitCell,opt);
pdf = @(e) exp(-beta * e);

%% initialising...
energy = zeros(nsamples, 1);
acceptance = zeros(nsamples, 1);
smpl = zeros(nsamples, size(T,2));
deformation = zeros(nsamples, 1);
deformationV = zeros(nsamples, size(T,1));
maxTheta = zeros(nsamples, 1);

%% main loop
u = u0(:);
for ct = 1:nsamples
    [e,~,~,~,theta] = Energy(u, T, extrudedUnitCell, opt);
    candidate = u + randn(size(T,2), 1) * delta;
    dTheta = theta - extrudedUnitCell.theta;
    % make sure that no faces are crossing over each other
    while ~isempty(find(abs(theta)>pi, 1))
        candidate = u + randn(size(T,2), 1) * delta;
    end
    
    [can_e,~,~,~,can_theta] = Energy(candidate, T, extrudedUnitCell, opt);
    acceptance(ct) = min(1, pdf(can_e)/pdf(e));
    
    % accept new candidate
    if rand < acceptance(ct)
        u = candidate;
        e = can_e;
        dTheta = can_theta - extrudedUnitCell.theta;
    end
    smpl(ct,:) = u';
    energy(ct) = e;
    V = T * u;
    deformationV(ct,:) = V;
    deformation(ct) = norm(V);
    [~,maxIdx] = max(abs(dTheta));
    maxTheta(ct) = dTheta(maxIdx);
end

if strcmp(opt.saveCSV, 'on')
    header = strcat('mh() beta=', num2str(opt.beta), ...
        ' delta=', num2str(opt.delta), ' Kface=', num2str(opt.Kface), ...
        ' Khinge=', num2str(opt.Khinge), ' Kedge=', num2str(opt.Kedge));
    saveCSV(opt, header, 'mh', 'smpl', smpl, 'acceptance', acceptance,...
        'energy', energy, 'deformationV', deformationV, 'maxTheta', maxTheta);
end


% %% mh & mh_plot example
% opt=initOpt('inputType','individual',...
%             'template','truncated tetrahedron',...
%             'plot','result',... %'scale', 1,... % deformation
%             'saveMovie', 'on', 'safeMovieAntiAlias', 0, ...
%             'saveCSV', 'on', 'saveGraph', 'on', ... graph from mh_plot
%             'beta', 1, 'delta', .02, ...
%             'interval',1,'saveFig','off','periodic','off',... 
%             'Khinge',0,'Kedge',1,'Kface',100,'KtargetAngle',1,...
%             'date', datestr(now, 'mmm-dd-yyyy'),...
%             'time', datestr(now,'HH-MM-SS'));
% [unitCell,extrudedUnitCell,opt]=buildGeometry(opt);
% [T]=linearConstr(unitCell,extrudedUnitCell,opt);
% u0 = zeros(1, size(T,2)); % starting point
% nsamples = uint32(2e2);
% tic;
% [smpl, acceptance, energy, deformationV, maxTheta] = mh(opt, u0, nsamples);
% toc;
% max(maxTheta)
% tic;
% mh_plot(nsamples, energy, deformationV, opt, acceptance, maxTheta)
% toc;