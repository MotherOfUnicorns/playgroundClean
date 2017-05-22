function [energy, deformationV] = grad_descent_num(opt, u0, max_iter)
% gradient descent using analytical jacobian
% u0 contains the coordinates of only the master nodes


[unitCell,extrudedUnitCell,opt]=buildGeometry(opt);
[T]=linearConstr(unitCell,extrudedUnitCell,opt);

eps = .0001;
func = @(u) Energy(u(:), T, extrudedUnitCell, opt);
helper = @(t,u) func(u);
jac_helper = @(t,u,e) numjac(helper, 0, u, e, eps); % u, e are the
        % deformation and energy at which the jacobian is evaluated
jac_func = @(u,e) jac_helper(0, u, e);


% variables to store
energy = zeros(max_iter, 1);
deformationV = zeros(max_iter, size(T,1));


% getting ready for the loop
tol = .1;
step = .0001;
u = u0(:);
ct = 1;
difference = 100;
while ct <= max_iter && difference > tol
    % energy and deformation of current u
    [e,~,~,~,theta] = func(u);
    energy(ct) = e;
    V = T * u;
    deformationV(ct, :) = V;

    % get the new solution from jacobian
    u_new = u - step * jac_func(u,e)';
    
    % calculate new deformation and 
    % difference between two consecutive solutions
    V_new = T * u_new;
    difference = norm(V - V_new); % ??? use this or simply norm(u - u_new)?
    
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
    %deformation = deformation(1:ct, :);
    deformationV = deformationV(1:ct, :);
end

% save data in csv file
if strcmp(opt.saveCSV, 'on')
    header = strcat('grad_descent_num() beta=', num2str(opt.beta), ...
        ' delta=', num2str(opt.delta), ' Kface=', num2str(opt.Kface), ...
        ' Khinge=', num2str(opt.Khinge), ' Kedge=', num2str(opt.Kedge));
    saveCSV(opt, header, 'grad_descent_num', 'energy', energy,...
        'deformationV', deformationV);
end


%% example grad_descent_num()
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
% max_iter = 1e2-1;
% hold off;
% 
% u0 = randn(size(T,2), 1) .* .5;
% [energy, deformationV] = grad_descent_num(opt, u0, max_iter);
% mh_plot(length(energy), energy, deformationV, opt)