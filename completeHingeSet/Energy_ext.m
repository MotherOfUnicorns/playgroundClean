function [E, dE,Eedge,Eface,Ehinge,theta]=Energy_ext(u,extrudedUnitCell,opt)
% extended energy function that also returns the analytical jacobian

dE=zeros(3*size(extrudedUnitCell.node,1),1);
Eedge=0;
Eface=0;

%APPLY DEFORMATION NODES
extrudedUnitCell.node=extrudedUnitCell.node+[u(1:3:end) u(2:3:end) u(3:3:end)];

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


%TOTAL ENERGY
E=Eedge+Eface+Ehinge;


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
    Jedge(i,3*extrudedUnitCell.edge(i,1)-2:3*extrudedUnitCell.edge(i,1))=-dx/L;
    Jedge(i,3*extrudedUnitCell.edge(i,2)-2:3*extrudedUnitCell.edge(i,2))=dx/L;
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
        Jface(rep,3*extrudedUnitCell.face{i}(1)-2:3*extrudedUnitCell.face{i}(1))=cross((coor3-coor2),(coor3-coor4));
        Jface(rep,3*extrudedUnitCell.face{i}(2)-2:3*extrudedUnitCell.face{i}(2))=cross((coor3-coor1),(coor4-coor1));
        Jface(rep,3*extrudedUnitCell.face{i}(3)-2:3*extrudedUnitCell.face{i}(3))=cross((coor4-coor1),(coor2-coor1));
        Jface(rep,3*extrudedUnitCell.face{i}(3+j)-2:3*extrudedUnitCell.face{i}(3+j))=cross((coor2-coor1),(coor3-coor1));
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
Jhinge=zeros(size(extrudedUnitCell.nodeHingeEx,1),size(extrudedUnitCell.node,1)*3);
for i=1:size(extrudedUnitCell.nodeHingeEx,1)
    extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,:),:);
    index(1:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:)-2;
    index(2:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:)-1;
    index(3:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:);
    [Jhinge(i,index),theta(i)]=JacobianHinge(extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,:),:));
end