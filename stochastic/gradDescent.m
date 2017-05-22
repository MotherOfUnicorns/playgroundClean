function [energy, deformationV] = gradDescent(opt, u0, max_iter)
% [energy, deformationV] = gradDescent(opt, u0, max_iter)
% 
% Gradient descnet method used to find local minima in energy function
% starting from some configuration *u0*.
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% opt      - options
% u0       - starting point
% max_iter - maximum number of iterations
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% energy       - energy of unit cell at each step
% deformationV - matrix containing vectors of deformation in each sample
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% EXAMPLE
% beta = 10;
% delta = .01;
% kface = 100;
% opt=initOpt('inputType','individual',...
%             'template','truncated tetrahedron',...
%             'plot','result',... %'scale', 1,... % deformation
%             'saveMovie', 'on', 'safeMovieAntiAlias', 0, ...
%             'saveCSV', 'on', 'saveGraph', 'on', ... graph from mh_plot
%             'beta', beta, 'delta', delta, ...
%             'interval',1,'saveFig','off','periodic','off',... 
%             'Khinge',1,'Kedge',100,'Kface',kface,'KtargetAngle',1,...
%             'date', datestr(now, 'mmm-dd-yyyy'),...
%             'time', datestr(now,'HH-MM-SS'),...
%             'constrEdge', 'off', 'constrFace','off');
% [~,extrudedUnitCell,opt]=buildGeometry(opt);
% 
% max_iter = 1e2;
% hold off;
% u0 = randn(3*size(extrudedUnitCell.node,1), 1) .* .5;
% [energy, deformationV] = gradDescent(opt, u0, max_iter);
% mh_plot(length(energy), energy, deformationV, opt)
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on May 22, 2017
% yun


[~,extrudedUnitCell,opt]=buildGeometry(opt);

% initialising...
energy = zeros(max_iter, 1);
deformationV = zeros(max_iter, length(u0));


% getting ready for the loops
tol = .01; % tolerence; terminate when *difference >= tol*
step = .0001; % step size
u = u0(:);

% first step
ct = 1;
[e,~,~,~,~,~] = Energy(u, extrudedUnitCell, opt);
energy(ct) = e;
deformationV(ct, :) = u;

ct = 2;
difference = 1; % difference in energy between two consecutive iterations
while ct <= max_iter+1 && difference > tol
    % energy of current u
    [e, de,~,~,~,theta] = Energy(u, extrudedUnitCell, opt);
    energy(ct) = e;
    deformationV(ct, :) = u;

    % get the new solution from jacobian
    u_new = u - step .* de(:);
    
    % calculate energy difference between two consecutive solutions
    difference = abs(energy(ct) - energy(ct-1));
    % difference = norm(u_new - u);
    
    % stop loop if one of the hinge angles are over the allowed range
    if ~isempty(find(abs(theta)>pi, 1))
        break;
    end
    
    % prepare for next step
    ct = ct + 1;
    u = u_new;
end


% trim off the excess if convergence occurred befor max_iter
if ct < max_iter+1
    energy = energy(1:ct, :);
    deformationV = deformationV(1:ct, :);
end

% save data in csv file
if strcmp(opt.saveCSV, 'on')
    header = strcat('grad_descent() beta=', num2str(opt.beta), ...
        ' delta=', num2str(opt.delta), ' Kface=', num2str(opt.Kface), ...
        ' Khinge=', num2str(opt.Khinge), ' Kedge=', num2str(opt.Kedge));
    saveCSV(opt, header, 'grad_descent', 'energy', energy,...
        'deformationV', deformationV);
end


%% EXAMPLE: Metropolis-Hastings with gradient descent
% opt=initOpt('inputType','individual',...
%             'template','truncated tetrahedron',...
%             'plot','result',... %'scale', 1,... % deformation
%             'saveMovie', 'on', 'safeMovieAntiAlias', 0, ...
%             'saveCSV', 'on', 'saveGraph', 'on', ... graph from mh_plot
%             'beta', .5, 'delta', .01, ...
%             'interval',1,'saveFig','off','periodic','off',... 
%             'Khinge',-1,'Kedge',100,'Kface',100,'KtargetAngle',1,...
%             'date', datestr(now, 'mmm-dd-yyyy'),...
%             'time', datestr(now,'HH-MM-SS'),...
%             'constrEdge', 'off', 'constrFace','off');
% [~,extrudedUnitCell,opt]=buildGeometry(opt);
% u0 = zeros(1, 3*size(extrudedUnitCell.node,1)); % starting point
% nsamples = uint32(2e2);
% tic;
% [acceptance, energy, deformationV] = mh(opt, u0, nsamples);
% toc;
% tic;
% mh_plot(nsamples, energy, deformationV, opt, acceptance)
% toc;
% 
% %% gradient descent starts here
% close all
% %opt.Khinge = -1; % still use negative hinge
% opt.Khinge = 1; %use positive hinge
% opt.date = datestr(now, 'mmm-dd-yyyy');
% opt.time = datestr(now,'HH-MM-SS');
% u0 = deformationV(end,:);
% max_iter = 1e2;
% tic;
% [energy, deformationV] = gradDescent(opt, u0, max_iter);
% toc;
% tic;
% mh_plot(length(energy), energy, deformationV, opt);
% toc;