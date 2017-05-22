function [energy, deformationV] = grad_descent_fmincon(opt, u0)
% gradient descent using fmincon
% uses Energy_ext() to get both energy and its jacobian
% u0 contains the coordinates of all nodes
% Feb 14, 2017
% yun li

% initialising
[unitCell,extrudedUnitCell,opt]=buildGeometry(opt);
options=optimoptions('fmincon','GradConstr','on','GradObj','on',...
    'tolfun', 1e-5, 'tolx', 1e-5,...
    'tolcon', 1e-6, 'Display','off',...
    'DerivativeCheck','off','maxfunevals',100000,...
    'OutputFcn', @outfun);
    % can add ('Display', 'iter') to show intermediate steps
    % 'StepTolerance', 1e-5,
    % StepTolerance is a lower bound on the size of a step
    % 'FunctionTolerance', 1e-3,
    % FunctionTolerance is a lower bound on the change in the value
    % of the objective function during a step
    
    
% Set up shared variables with OUTFUN
history.x = [];
history.fval = [];

% define constrains
[Aeq, Beq]=linearConstraint(unitCell,extrudedUnitCell,opt);

% define problem with problem struct
problem.options = options;
problem.solver = 'fmincon';
problem.objective = @(u)Energy_ext(u, extrudedUnitCell,opt);
problem.x0 = u0;
problem.Aineq = [];
problem.bineq = [];
problem.Aeq = Aeq;
problem.beq = Beq;
problem.lb = [];
problem.ub = [];
problem.nonlcon = @(u)nonlinearConstr(u,extrudedUnitCell,opt);
%problem.nonlcon = [];

% find local minima with fmincon
fmincon(problem);
energy = history.fval;
deformationV = reshape(history.x, length(energy),[]);

% save file as csv
if strcmp(opt.saveCSV, 'on')
    header = strcat('grad_descent_fmincon() beta=', num2str(opt.beta), ...
        ' delta=', num2str(opt.delta), ' Kface=', num2str(opt.Kface), ...
        ' Khinge=', num2str(opt.Khinge), ' Kedge=', num2str(opt.Kedge));
    saveCSV(opt, header, 'grad_descent_fmincon', 'energy', energy,...
        'deformationV', deformationV);
end


%%%%%%%%%%%%%%%%%%%%%%%%%
    function stop = outfun(x,optimValues,state)
    % output function for fmincon

    stop = false;
    switch state
        case 'init'
            hold on
        case 'iter'
            % Concatenate current point and objective function
            % value with history. x must be a row vector.
            history.fval = [history.fval; optimValues.fval];
            history.x = [history.x; x];
        case 'done'
            hold off
        otherwise
    end
    end

end


function [Aeq, Beq]=linearConstraint(unitCell,extrudedUnitCell,opt)
% linear constrains

%FIX NODE CONSTRAINTS
%IMPROVE FOLLOWING - AUTOMATIC DEPENDING ON NODES OF FACE 1
nodeFix=extrudedUnitCell.face{1};
e1=extrudedUnitCell.node(nodeFix(2),:)-extrudedUnitCell.node(nodeFix(1),:);
e2=extrudedUnitCell.node(nodeFix(3),:)-extrudedUnitCell.node(nodeFix(2),:);
e1=e1/norm(e1);
e2=e2/norm(e2);
e3=cross(e2,e1);
e3=e3/norm(e3);

Aeq=zeros(6,3*size(extrudedUnitCell.node,1));
Beq=zeros(6,1);
Aeq(1,3*nodeFix(2)-2)=1;
Aeq(2,3*nodeFix(2)-1)=1;
Aeq(3,3*nodeFix(2))=1;
Aeq(4,3*nodeFix(1)-2:3*nodeFix(1))=e3;
Aeq(5,3*nodeFix(1)-2:3*nodeFix(1))=e2;
Aeq(6,3*nodeFix(3)-2:3*nodeFix(3))=e3;

%MERGE NODES AT INITIALLY SAME LOCATION
rep=size(Aeq,1);
for i=1:size(unitCell.internalFacePairs,1)   
    for j=1:length(unitCell.Polyhedron(unitCell.internalFacePairs(i,2)).faceNodeExtrude{unitCell.internalFacePairs(i,1)})
        index1=unitCell.Polyhedron(unitCell.internalFacePairs(i,2)).faceNodeExtrude{unitCell.internalFacePairs(i,1)}(j);
        for k=1:length(unitCell.Polyhedron(unitCell.internalFacePairs(i,4)).faceNodeExtrude{unitCell.internalFacePairs(i,3)})
            index2=unitCell.Polyhedron(unitCell.internalFacePairs(i,4)).faceNodeExtrude{unitCell.internalFacePairs(i,3)}(k);
            if norm(extrudedUnitCell.node(index2,:)'-extrudedUnitCell.node(index1,:)')<opt.Lextrude/1e6
                rep=rep+1;
                Aeq(3*rep-2:3*rep,:)=zeros(3,size(extrudedUnitCell.node,1)*3);
                Aeq(3*rep-2:3*rep,3*index1-2:3*index1)=[1 0 0; 0 1 0; 0 0 1];
                Aeq(3*rep-2:3*rep,3*index2-2:3*index2)=[-1 0 0; 0 -1 0; 0 0 -1];
                Beq(3*rep-2:3*rep)=0;
            end
        end
     end
end

%PERIODIC NODAL CONSTRAINTS
if strcmp(opt.periodic,'on')
    nref=length(extrudedUnitCell.ref);
    rep=size(Aeq,1);
    for i=1:size(unitCell.possibleAlpha,1)
        for j=1:size(extrudedUnitCell.node,1)-nref
            coor1=extrudedUnitCell.node(j,:)';
            for k=1:size(extrudedUnitCell.node,1)-nref
                coor2=extrudedUnitCell.node(k,:)';
                if norm(coor2-coor1-unitCell.l'*unitCell.possibleAlpha(i,:)')<1e-6
                    rep=rep+1;
                    %sprintf('%d, node 1 = %d, node 2 =%d',[rep,j,k])
                    Aeq(3*rep-2:3*rep,:)=zeros(3,size(extrudedUnitCell.node,1)*3);
                    Aeq(3*rep-2:3*rep,3*j-2:3*j)=[1 0 0; 0 1 0; 0 0 1];
                    Aeq(3*rep-2:3*rep,3*k-2:3*k)=[-1 0 0; 0 -1 0; 0 0 -1];
                    for l=1:nref
                        Aeq(3*rep-2:3*rep,3*extrudedUnitCell.ref(l)-2:3*extrudedUnitCell.ref(l))=unitCell.possibleAlpha(i,l)*[-1 0 0; 0 -1 0; 0 0 -1];
                    end
                    Beq(3*rep-2:3*rep)=0;
                end
            end
        end
    end
end
end

function [C,Ceq,DC,DCeq]=nonlinearConstr(u,extrudedUnitCell,opt)
% nonlinear constraints
C1=[]; 
DC1=[]; 
Ceq1=[]; Ceq2=[];
DCeq1=[]; DCeq2=[];

%APPLY DEFORMATION NODES
extrudedUnitCell.node=extrudedUnitCell.node+[u(1:3:end) u(2:3:end) u(3:3:end)];

%%%%%%%%%%%%%%%%%%%%%%
%INEQUALITY CONSTRAINS
%%%%%%%%%%%%%%%%%%%%%%
%MAXIMUM AND MINIMUM ANGLES
[C, DC]=getHinge(extrudedUnitCell);
C=[-C-pi*opt.constAnglePerc; C-pi*opt.constAnglePerc];
DC=[-DC; DC]';

%%%%%%%%%%%%%%%%%%%%%%
%INEQUALITY CONSTRAINS
%%%%%%%%%%%%%%%%%%%%%%
%CONSTRAINT ASSOCIATED TO EDGE STRETCHING
if strcmp(opt.constrEdge,'on')
    [Ceq1, DCeq1]=getEdge(extrudedUnitCell);
end
%ENERGY ASSOCIATED TO FACE BENDING
if strcmp(opt.constrFace,'on')
    [Ceq2, DCeq2]=getFace(extrudedUnitCell);
end
%FINAL CONSTRAINTS
Ceq=[Ceq1; Ceq2];
DCeq=[DCeq1; DCeq2]';
end


%% example grad_descent_fmincon
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
%             'constrEdge', 'off', 'constrFace','off',...
%             'constAnglePerc',0.98);
% [unitCell,extrudedUnitCell,opt]=buildGeometry(opt);
% [T]=linearConstr(unitCell,extrudedUnitCell,opt);
% u0 = randn(size(T,1), 1) .* .1;
% tic;
% [energy, deformationV] = grad_descent_fmincon(opt,u0);
% toc;
% mh_plot(length(energy), energy, deformationV, opt)