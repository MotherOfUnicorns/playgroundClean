function [smpl, accept] = mh_plot(beta, delta, nupdate, nsamples, varargin)
% simulates with Metropolis-Hastings algorithm and also plots the results

% initialising...
opt=initOpt('inputType','individual',...
    'template','truncated tetrahedron',...
    'plot','result',...
    'interval',1,'saveFig','off','periodic','off',...
    'Khinge',0.0005,'Kedge',1,'Kface',1,'KtargetAngle',1);
[unitCell,extrudedUnitCell,opt]=buildGeometry(opt);
[T]=linearConstr(unitCell,extrudedUnitCell,opt);

u0 = zeros(1, size(T,2));
pdf = @(u) exp(-beta * Energy(u(:),T,extrudedUnitCell,opt));
smpl = zeros(nsamples + 1, size(u0,2)); % include inital condition in results
accept = zeros(nsamples + 1, 1);


% initialising stuff to display in plots
ct = 0;
smpl(1, :) = u0(:);
energy = zeros(floor(nsamples/nupdate), 1);
energy(0) = Energy(u0(:),T,extrudedUnitCell,opt);
deformation = zeros(floor(nsamples/nupdate), 1);
deformation(0) = T * u0(:);


% initialising plots
figure; hold on;
subplot(3,2,1);
h1 = plot(ct, accept(ct + 1));
xlabel('n')
ylabel('acceptance')

subplot(3,2,3);
h2 = plot(ct, energy(ct + 1));
xlabel('n')
ylabel('energy')

subplot(3,2,4);
h3 = histpgram(energy, floor(nsamples/100));
set(gca,'YDir','reverse')
view(-90,90) % swap x and y axis
xlabel('count')
ylabel('energy')

subplot(3,2,5);
h4 = plot(ct, deformation(ct + 1));
xlabel('n')
ylabel('deformation')

subplot(3,2,6);
h5 = plot(energy, deformation);
xlabel('energy')
ylabel('deformation')

figure; hold on;
result.numMode = 1;
outputResults(unitCell,extrudedUnitCell,result,opt);
hCell = gcf; % do something here


u = u0(:);
for ct = 1:nsamples
    
    candidate = u + randn(size(u, 1), 1) * delta;
    acceptance = min(1, pdf(candidate)/pdf(u));
    if rand < acceptance
        u = candidate;
    end
    smpl(ct + 1,:) = u';
    accept(ct + 1) = acceptance;
    
    
    % update figures every n steps
    if mod(ct, nupdate) == 0
        figIdx = int(ct/update) + 1;% index for updating figures
        
        
        % update view of unit cell
        V = T * u'; % deformation
        result.numMode=1;
        energy(figIdx) = Energy(u', T, extrudedUnitCell, opt);
        result.deform.V=[V(1:3:end) V(2:3:end) V(3:3:end)];
        result.deform.Ve=V(:);
        
        deformation(figIdx) = norm(V);
        %outputResults(unitCell,extrudedUnitCell,result,opt)
        
       % update plots
        set(h1, 'ydata', accept)
        set(h2, 'xdata', 1:nupdate:ct, 'ydata', energy(1:figIdx))
        set(h3) % update histogram?
        set(h4, 'xdata', 1:nupdate:ct, 'ydata', deformation(1:figIdx))
        set(h5, 'xdata', energy, 'ydata', energy)
        set(hCell)
        drawnow;
    end
    
    
end


