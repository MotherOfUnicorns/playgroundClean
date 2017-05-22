% all functions (mh.m, annealing.m, breathing.m, gradDescent.m) developed
% for STOCHASTIC are presented here.
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on May 22, 2017
% yun


close all; clear

opt=initOpt('inputType','individual',...
        'template','truncated tetrahedron',...
        'plot','result',... %'scale', 1,... % deformation
        'saveMovie', 'on', 'safeMovieAntiAlias', 0, ...
        'saveCSV', 'on', 'saveGraph', 'on', ... graph from mh_plot
        'beta', .1, 'delta', .005, ...
        'interval',1,'saveFig','off','periodic','off',... 
        'Khinge',-1,'Kedge',100,'Kface',100,'KtargetAngle',1,...
        'date', datestr(now, 'mmm-dd-yyyy'),...
        'time', datestr(now,'HH-MM-SS'),...
        'constrEdge', 'off', 'constrFace','off',...
        'constAnglePerc',0.99);
[~,extrudedUnitCell,opt]=buildGeometry(opt);
u0 = zeros(1, 3*size(extrudedUnitCell.node,1)); % starting point
nsamples = 1e2;
nrepeats = 2; % only useful for breathing modes

func = 'breathing_gradient';
% *func* must be one of the following:
% mh                 - Metropolis-Hastings 
% annealing          - simulated annealing
% breathing          - Metropolis-Hastings with alternating high/low T
% breathing_gradient - breathing mode plus optimisation by gradDescent


switch func
    case 'mh'
        [acceptance,energy,deformationV,maxTheta] =  mh(opt,u0,nsamples);
        mh_plot(nsamples, energy, deformationV, opt, acceptance, maxTheta)
    
    case 'annealing'
        [acceptance, energy, deformationV, maxTheta] = ...
            annealing(opt, u0, nsamples);
        mh_plot(nsamples, energy, deformationV, opt, acceptance, maxTheta)
    
    case 'breathing'
        [acceptance, energy, deformationV, maxTheta] = ...
            breathing(opt, u0, nsamples, nrepeats);
        mh_plot(nsamples*nrepeats, energy, deformationV, opt, ...
            acceptance, maxTheta)
    
    case 'breathing_gradient'
        [acceptance1, energy1, deformationV1, maxTheta1] = ...
            breathing(opt, u0, nsamples, nrepeats);
        mh_plot(nsamples*nrepeats, energy1, deformationV1, opt, ...
            acceptance1, maxTheta1)
        
        % start optimsations
        [~,idx] = max(maxTheta1);
        u0 = deformationV1(idx, :)';
        opt.Khinge = 1;
        
        
        % feed result to gradient descent for further optimisation
        max_iter = 1e4;
        opt.date = datestr(now, 'mmm-dd-yyyy');
        opt.time = datestr(now,'HH-MM-SS');
        [energy2, deformationV2] = gradDescent(opt,u0,max_iter);
        close all;
        mh_plot(length(energy2), energy2, deformationV2, opt);
    otherwise
        error('Function not defined')
end