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