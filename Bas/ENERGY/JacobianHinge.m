function [t]=JacobianHinge(p0)

%Node coordinates
p=reshape(p0',12,1);
%Vectors
a=p(4:6)-p(1:3);    %Rotation axis
b=p(7:9)-p(1:3);    %Vector on face 1
c=p(10:12)-p(1:3);  %Vector on face 2
 
ka=norm(a);
kab=sqrt((a'*a)*(b'*b)-(a'*b)^2);
kca=sqrt((c'*c)*(a'*a)-(c'*a)^2);
na=a/ka;

kab=real(kab);
kca=real(kca);
nab=crossvector(a,b)/kab;
nca=crossvector(c,a)/kca;

detf=na'*(crossvector(nab,nca));
dotf=nab'*nca;
t = atan2(detf,dotf);

function w=crossvector(u,v)
w=[u(2,:).*v(3,:)-u(3,:).*v(2,:);-u(1,:).*v(3,:)+u(3,:).*v(1,:);u(1,:)*v(2,:)-u(2,:)*v(1,:)];