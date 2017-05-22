function [result,extrudedUnitCell, opt] = ...
    findDeformation(unitCell,extrudedUnitCell,opt)
% first with fmincon, then with gradDescent (with penalty for faces
% crossing over)
% last modified on Mar 31, 2017

% save this for second optimisation with gradient descent
opt.extrudedUnitCell = extrudedUnitCell;

%Show details geometries (if requested)
if strcmp(opt.plot,'info')
    result=[];
else
	% first use fmincon
    [V1, extrudedUnitCell, exitFlag] = ...
        nonlinearFolding(unitCell, extrudedUnitCell, opt);
    opt.exitFlag1 = exitFlag;
    
    % then use gradDescent
    max_iter = 100000;
    u0 = V1(:,end);
    [V2, exitFlag] = gradDescent(opt, u0, max_iter);
    opt.exitFlag2 = exitFlag;
    
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
%ENERGY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E, dE,Eedge,Eface,Ehinge,EtargetAngle,theta] = ...
    Energy(u,extrudedUnitCell,opt)

E=0; dE=zeros(3*size(extrudedUnitCell.node,1),1);
Eedge=0;
Eface=0;
Ehinge=0;
EtargetAngle=0;

%APPLY DEFORMATION NODES
extrudedUnitCell.node = extrudedUnitCell.node ...
                        +[u(1:3:end) u(2:3:end) u(3:3:end)];

%ENERGY ASSOCIATED TO EDGE STRETCHING
if strcmp(opt.constrEdge,'off')
    [dEdge, Jedge]=getEdge(extrudedUnitCell);
    Eedge=1/2*opt.Kedge*sum(dEdge.^2);
    dE=dE+opt.Kedge*Jedge'*dEdge;
end

%ENERGY ASSOCIATED TO FACE BENDING
if strcmp(opt.constrFace,'off')
    [dFace, Jface]=getFace(extrudedUnitCell);
    Eface=opt.Kface/2*sum(dFace.^2);
    dE=dE+opt.Kface*Jface'*dFace;
end

%ENERGY ASSOCIATED TO HINGE BENDING
[theta, Jhinge]=getHinge(extrudedUnitCell);
Ehinge=1/2*opt.Khinge*sum((theta-extrudedUnitCell.theta).^2);
dE=dE+opt.Khinge*(Jhinge'*(theta-extrudedUnitCell.theta));

%ENERGY ASSOCIATED TO TARGET HINGE ANGLES
if size(extrudedUnitCell.angleConstr,1)==0
    dtheta=[];
else
    dtheta = theta(extrudedUnitCell.angleConstr(:,1))...
            -extrudedUnitCell.angleConstr(:,2);
    dE = dE + opt.KtargetAngle...
        *Jhinge(extrudedUnitCell.angleConstr(:,1),:)' * dtheta;
end
EtargetAngle=1/2*opt.KtargetAngle*sum(dtheta.^2);

%TOTAL ENERGY
E=Eedge+Eface+Ehinge+EtargetAngle;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EDGE LENGTH AND JACOBIAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dEdge, Jedge]=getEdge(extrudedUnitCell)
dEdge=zeros(size(extrudedUnitCell.edge,1),1);
Jedge=zeros(length(extrudedUnitCell.edge),size(extrudedUnitCell.node,1)*3);   
for i=1:size(extrudedUnitCell.edge,1)
    coor1=extrudedUnitCell.node(extrudedUnitCell.edge(i,1),:);
    coor2=extrudedUnitCell.node(extrudedUnitCell.edge(i,2),:);
    dx=coor2-coor1;
    L=sqrt(dx*dx');
    dEdge(i)=L-extrudedUnitCell.edgeL(i);            
    Jedge(i,3*extrudedUnitCell.edge(i,1)-2:3*extrudedUnitCell.edge(i,1))...
        = -dx/L;
    Jedge(i,3*extrudedUnitCell.edge(i,2)-2:3*extrudedUnitCell.edge(i,2))...
        = dx/L;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FACE OUT OF PLANE DEFORMATION AND JACOBIAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dFace, Jface]=getFace(extrudedUnitCell)
rep=0;
tnf=0;
for i=1:length(extrudedUnitCell.face)
    tnf=tnf+length(extrudedUnitCell.face{i})-3;
end
dFace=zeros(tnf,1);
Jface=zeros(tnf,size(extrudedUnitCell.node,1)*3);
for i=1:length(extrudedUnitCell.face)
    coor1=extrudedUnitCell.node(extrudedUnitCell.face{i}(1),:);
    coor2=extrudedUnitCell.node(extrudedUnitCell.face{i}(2),:);
    coor3=extrudedUnitCell.node(extrudedUnitCell.face{i}(3),:);
    a=cross(coor2-coor1,coor3-coor1);
    for j=1:length(extrudedUnitCell.face{i})-3
        rep=rep+1;
        coor4=extrudedUnitCell.node(extrudedUnitCell.face{i}(3+j),:);
        Jface(rep,3*extrudedUnitCell.face{i}(1)-2:...
            3*extrudedUnitCell.face{i}(1))...
            = cross((coor3-coor2),(coor3-coor4));
        Jface(rep,3*extrudedUnitCell.face{i}(2)-2:...
            3*extrudedUnitCell.face{i}(2))...
            =cross((coor3-coor1),(coor4-coor1));
        Jface(rep,3*extrudedUnitCell.face{i}(3)-2:...
            3*extrudedUnitCell.face{i}(3)) = ...
            cross((coor4-coor1),(coor2-coor1));
        Jface(rep,3*extrudedUnitCell.face{i}(3+j)-2:...
            3*extrudedUnitCell.face{i}(3+j)) ...
            =cross((coor2-coor1),(coor3-coor1));
        dFace(rep)=(coor4-coor1)*a';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HINGE ANGLE AND JACOBIAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta, Jhinge]=getHinge(extrudedUnitCell)
theta=zeros(size(extrudedUnitCell.nodeHingeEx,1),1);
Jhinge = zeros(size(extrudedUnitCell.nodeHingeEx,1),...
               size(extrudedUnitCell.node,1)*3);
for i=1:size(extrudedUnitCell.nodeHingeEx,1)
    extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,:),:);
    index(1:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:)-2;
    index(2:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:)-1;
    index(3:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:);
    [Jhinge(i,index),theta(i)] = ...
        JacobianHinge(extrudedUnitCell.node(...
                    extrudedUnitCell.nodeHingeEx(i,:),:));
end




function [deformationV, exitFlag] = gradDescent(opt, u0, max_iter)
% gradient descent method without using fmincon, a penalty is imposed for
% preventing faces corssing over each other.
% u0 contains the coordinates of all nodes
% Apr 04, 2017
% yun li

%unitCell = opt.unitCell;
extrudedUnitCell = opt.extrudedUnitCell;
extrudedUnitCell.angleConstr = [];
exitFlag = 1;

% getting ready for the loop
tol = 1e-9;
step = 1e-1;
u = u0(:);
ct = 0;
difference = Inf;
e_prev = Inf;
while ct <= max_iter && difference > tol
    % energy of current u
    [e, de, ~,~,~,~, theta] = Energy(u, extrudedUnitCell, opt);
    % uses penalty to prevent faces from crossing over
    e = e + penalty(theta);
    
%     % stop loop if one of the hinge angles are over the allowed range
%     if ~isempty(find(abs(theta)>pi, 1))
%         warning('faces overlap!')
%         break;
%     end
    
    % get the new solution from jacobian
    u_new = u - step .* de(:);
    
    % calculate difference between two consecutive solutions
    difference = abs(e - e_prev) / ...
        size(extrudedUnitCell.edge,1); % normalised to the number of hinges
    
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
% a penalty function which allows values of theta that are in [-pi, pi]
ePen = sum(max(0, 0.5 .* (theta - pi) .* (theta + pi)));