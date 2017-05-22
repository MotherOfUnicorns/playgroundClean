function [T]=linearConstr(unitCell,extrudedUnitCell,opt)

%FIX NODE CONSTRAINTS
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

%DETERMINE TRANSFORMATION MATRIX BASED ON CONSTRAINTS
[A,indexSlave]=frref(sparse(Aeq),'1e-8','s'); %Reduced row echelon form
indexMaster=setdiff(1:size(A,2),indexSlave);
indexConstr=1:length(indexSlave);
Am=A(indexConstr,indexMaster);
T=diag(ones(3*size(extrudedUnitCell.node,1),1));
T(indexSlave,indexMaster)=T(indexSlave,indexMaster)-Am;
T=T(:,indexMaster);