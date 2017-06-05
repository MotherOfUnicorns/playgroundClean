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






