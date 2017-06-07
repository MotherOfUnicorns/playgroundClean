function [result, extrudedUnitCell, exitFlag1, exitFlag2, opt] = ...
    findDeformation(unitCell, extrudedUnitCell, opt)
% Two steps of optimisation: first with fmincon, then with gradDescent
% (with penalty for faces crossing over)
% Adapted from Bas' code
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Jun 05, 2017
% yun

% save this for second optimisation with gradient descent
original_extrudedUnitCell = extrudedUnitCell;

%Show details geometries (if requested)
if strcmp(opt.plot,'info')
    exitFlag1 = 0;
    exitFlag2 = 0;
    result=[];
else % optimisations...
	% first use fmincon
    [V1, extrudedUnitCell, exitFlag1] = ...
        nonlinearFolding(unitCell, extrudedUnitCell, opt);

    % then use gradDescent
    max_iter = 100000;
    u0 = V1(:,end);
    [V2, exitFlag2] = ...
        gradDescent(original_extrudedUnitCell, opt, u0, max_iter);
    
    % combine results
    V = [V1, V2];
    result.numMode=size(V,2);

    %Create different way of storing results
    for ct=1:result.numMode
        result.deform(ct).V = [V(1:3:end,ct) V(2:3:end,ct) V(3:3:end,ct)];
        result.deform(ct).Ve = V(:,ct);
        %result.deform(ct).interV(1).V = result.deform(ct).V;
        %result.deform(ct).interV(1).Ve=V(:,ct);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NON-LINEAR ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V,extrudedUnitCell, exitFlag] =  ...
    nonlinearFolding(unitCell, extrudedUnitCell, opt)

%INITIALIZE LINEAR CONSTRAINTS and output
[Aeq, Beq]=linearConstr(unitCell,extrudedUnitCell,opt);
V = zeros(3*size(extrudedUnitCell.node,1), 2);
u0 = V(:, 1);
theta0 = extrudedUnitCell.theta;
extrudedUnitCell.angleConstr = opt.angleConstrFinal(1).val;

% here starts the main thing
tic
% if size(opt.angleConstrFinal(1).val,1)~=0
%     extrudedUnitCell.angleConstr
%     
%     
%     extrudedUnitCell.angleConstr(:,2) = ...
%         theta0(extrudedUnitCell.angleConstr(:,1))...
%         -theta0(extrudedUnitCell.angleConstr(:,1)...
%         +opt.angleConstrFinal(1).val(:,2));
% end
%Determine new equilibrium
[V(:,2),~,exitFlag]=fmincon(@(u) ...
    Energy(u,extrudedUnitCell,opt),u0,[],[],...
    Aeq,Beq,[],[],@(u) ...
    nonlinearConstr(u,extrudedUnitCell,opt),opt.options);

%u0=V(:,2);
%Determine energy of that equilibrium
% [E(2,iter),~,Eedge(2,iter),Eface(2,iter),...
%     Ehinge(2,iter),EtargetAngle(2,iter)] = ...
%     Energy(u0,extrudedUnitCell,opt);
% t=toc;
% fprintf(', time: %1.2f, exitflag: %d\n',t,exfl)
%end
%update angle starting point
%[~,~,~,~,~,~,theta0]=Energy(u0,extrudedUnitCell,opt);
%end

% result.E=E;
% result.Eedge=Eedge;
% result.Eface=Eface;
% result.Ehinge=Ehinge;
% result.EtargetAngle=EtargetAngle;
% result.Theta = theta0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NONLINEAR CONSTRAINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C,Ceq,DC,DCeq]=nonlinearConstr(u,extrudedUnitCell,opt)

C1=[]; 
DC1=[]; 
Ceq1=[]; Ceq2=[];
DCeq1=[]; DCeq2=[];

%APPLY DEFORMATION NODES
extrudedUnitCell.node = extrudedUnitCell.node ...
                        +[u(1:3:end) u(2:3:end) u(3:3:end)];

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LINEAR CONSTRAINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Aeq, Beq]=linearConstr(unitCell,extrudedUnitCell,opt)

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
    for j=1:length(unitCell.Polyhedron(unitCell.internalFacePairs(i,2))...
            .faceNodeExtrude{unitCell.internalFacePairs(i,1)})
        index1 = unitCell.Polyhedron(unitCell.internalFacePairs(i,2))...
            .faceNodeExtrude{unitCell.internalFacePairs(i,1)}(j);
        for k=1:length(unitCell.Polyhedron(unitCell...
                .internalFacePairs(i,4)).faceNodeExtrude{unitCell...
                .internalFacePairs(i,3)})
            index2=unitCell.Polyhedron(unitCell.internalFacePairs(i,4))...
                .faceNodeExtrude{unitCell.internalFacePairs(i,3)}(k);
            if norm(extrudedUnitCell.node(index2,:)' ...
                    -extrudedUnitCell.node(index1,:)') < opt.Lextrude/1e6
                rep=rep+1;
                Aeq(3*rep-2:3*rep,:) = ...
                    zeros(3,size(extrudedUnitCell.node,1)*3);
                Aeq(3*rep-2:3*rep,3*index1-2:3*index1) = ...
                    [1 0 0; 0 1 0; 0 0 1];
                Aeq(3*rep-2:3*rep,3*index2-2:3*index2) = ...
                    [-1 0 0; 0 -1 0; 0 0 -1];
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
                if norm(coor2-coor1-unitCell.l' ...
                        *unitCell.possibleAlpha(i,:)') < 1e-6
                    rep=rep+1;
                    %sprintf('%d, node 1 = %d, node 2 =%d',[rep,j,k])
                    Aeq(3*rep-2:3*rep,:) = ...
                        zeros(3,size(extrudedUnitCell.node,1)*3);
                    Aeq(3*rep-2:3*rep,3*j-2:3*j) = ...
                        [1 0 0; 0 1 0; 0 0 1];
                    Aeq(3*rep-2:3*rep,3*k-2:3*k) = ...
                        [-1 0 0; 0 -1 0; 0 0 -1];
                    for l=1:nref
                        Aeq(3*rep-2:3*rep,3*extrudedUnitCell.ref(l)...
                            -2:3*extrudedUnitCell.ref(l)) = ...
                            unitCell.possibleAlpha(i,l) ...
                            *[-1 0 0; 0 -1 0; 0 0 -1];
                    end
                    Beq(3*rep-2:3*rep)=0;
                end
            end
        end
    end
end





function [deformationV, exitFlag] = ...
    gradDescent(original_extrudedUnitCell, opt, u0, max_iter)
% gradient descent method with a penalty imposed for preventing faces from
% corssing over each other.
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% original_extrudedUnitCell - extruded unit cell without applying any
%                             deformations
% opt      - options
% u0       - coordinates of all nodes
% max_iter - maximum number of iterations
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% deformationV - an array containing the deformation in each iteration
% exitFlag     - if exitFlag = 1, convergence occured within max_iter
%                if exitFlag = 0, max_iter is reached without convergence
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on May 22, 2017
% yun


original_extrudedUnitCell.angleConstr = [];
exitFlag = 1;

% getting ready for the loop
tol = 1e-9;
step = 1e-1; % step size
u = u0(:);
ct = 0;
difference = Inf;
e_prev = Inf;
while ct <= max_iter && difference > tol
    % energy of current u
    [e, de, ~,~,~,~, theta] = Energy(u, original_extrudedUnitCell, opt);
    % uses penalty to prevent faces from crossing over
    e = e + penalty(theta);
    
    % get the new solution from jacobian
    u_new = u - step .* de(:);
    
    % calculate difference between two consecutive solutions
    difference = abs(e - e_prev) / ...
        size(original_extrudedUnitCell.edge,1); % normalised to the number of hinges
    
    % prepare for next step
    ct = ct + 1;
    u = u_new;
    e_prev = e;
end

if ct > max_iter
    exitFlag = 0;
end

deformationV = u(:);


function ePen = penalty(theta)
% A penalty function which punished theta value outside the range of
% [-0.95*pi, 0.95*pi].
% The range is not [-pi, pi]
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on May 22, 2017
% yun

ePen = sum(max(0, 0.5 .* (theta - 0.95*pi) .* (theta + 0.95*pi)));