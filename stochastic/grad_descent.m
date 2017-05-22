function [energy, deformationV] = grad_descent(opt, u0, max_iter)
% naive gradient descent without using fmincon
% uses Energy_ext() to get both energy and its jacobian
% u0 contains the coordinates of all nodes
% Feb 10, 2017
% yun li


[unitCell,extrudedUnitCell,opt]=buildGeometry(opt);

% variables to store
energy = zeros(max_iter, 1);
deformationV = zeros(max_iter, length(u0));


% getting ready for the loop
tol = .01;
step = .0001;
u = u0(:);
ct = 1;
difference = 1;
while ct <= max_iter && difference > tol
    % energy of current u
    [e, de, ~,~,~, theta] = Energy_ext(u, extrudedUnitCell, opt);
    energy(ct) = e;
    deformationV(ct, :) = u;

    % get the new solution from jacobian
    u_new = u - step .* de(:);
    
    % calculate difference between two consecutive solutions
    difference = norm(u_new - u);
    
    % stop loop if one of the hinge angles are over the allowed range
    if ~isempty(find(abs(theta)>pi, 1))
        break;
    end
    
    % prepare for next step
    ct = ct + 1;
    u = u_new;
end


% trim off the excess
if ct < max_iter
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

%% example
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
% [unitCell,extrudedUnitCell,opt]=buildGeometry(opt);
% [T]=linearConstr(unitCell,extrudedUnitCell,opt);
% max_iter = 1e2;
% hold off;
% 
% u0 = randn(size(T,1), 1) .* .5;
% [energy, deformationV] = grad_descent(opt, u0, max_iter);
% mh_plot(length(energy), energy, deformationV, opt)


%% Metropolis-Hastings with gradient descent
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
% [unitCell,extrudedUnitCell,opt]=buildGeometry(opt);
% [T]=linearConstr(unitCell,extrudedUnitCell,opt);
% u0 = zeros(1, size(T,2)); % starting point
% nsamples = uint32(2e2);
% tic;
% [smpl, acceptance, energy, deformationV] = mh(opt, u0, nsamples);
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
% max_iter = 5e2;
% tic;
% [ energy, deformationV] = grad_descent(opt, u0, max_iter);
% toc;
% tic;
% mh_plot(length(energy), energy, deformationV, opt);
% toc;