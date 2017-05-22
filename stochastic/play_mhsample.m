
%%
opt=initOpt('inputType','individual',...
            'template','truncated tetrahedron',...
            'plot','result',...
            'saveMovie', 'on', 'safeMovieAntiAlias', 0, ...
            'interval',1,'saveFig','off','periodic','off',...
            'Khinge',0.0005,'Kedge',1,'Kface',1,'KtargetAngle',1);  
[unitCell,extrudedUnitCell,opt]=buildGeometry(opt);
[T]=linearConstr(unitCell,extrudedUnitCell,opt);

%%
%rng default  % For reproducibility
beta = 150; % some constant --> 1/temperature
delta = .01;
%u0 = rand(size(T,2), 1); % random starting point
u0 = zeros(1, size(T,2));
pdf = @(u) exp(-beta * Energy(u(:),T,extrudedUnitCell,opt));
proppdf = @(u,v) unifpdf(v-u, -delta, delta);
proprnd = @(u) u + (rand(1, size(T,2))*2 - ones(1, size(T,2))) * delta;
nsamples = 100;
[x, accept] = mhsample(u0, nsamples, 'pdf', pdf, 'proppdf', proppdf,...
   'proprnd',proprnd, 'symmetric', 1);


result.numMode=nsamples;
displacement = zeros(1,nsamples);
e = zeros(nsamples, 1);
for i = 1:nsamples
    V=T*x(i,:)'; % deformation
	e(i) = Energy(x(i,:)', T, extrudedUnitCell, opt);
    result.deform(i).V=[V(1:3:end) V(2:3:end) V(3:3:end)];
    result.deform(i).Ve=V(:);
    displacement(i) = norm(V);
end

figure;
title(strcat('beta', num2str(beta)))

subplot(3,2,1)
plot(accept, 'o')
xlabel('n')
ylabel('acceptance')

subplot(3,2,3)
plot(e, 'o');
xlabel('n')
ylabel('energy')

subplot(3,2,4)
histogram(e,1000)
set(gca,'YDir','reverse')
view(-90,90) % swap x and y axis
xlabel('energy')
ylabel('count')

subplot(3,2,5)
plot(displacement,'o')
xlabel('n')
ylabel('deformation')

subplot(3,2,6)
plot(displacement, e, 'o')
xlabel('energy')
ylabel('deformation')
%savefig(strcat('beta=', num2str(beta), 'nsample=', num2str(nsamples), '.fig'))

outputResults(unitCell,extrudedUnitCell,result,opt)
%savefig(strcat('beta=', num2str(beta), 'nsample=', num2str(nsamples), '_unitCell.fig'))







