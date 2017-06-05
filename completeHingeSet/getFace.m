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
