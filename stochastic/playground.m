% last modified on Feb 24, 2017
% yun
% everything aggregated together
% 
% FUNC(str) is one of the following functions:
%       mh; breathing; annealing; breathing_gradient

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
[unitCell,extrudedUnitCell,opt]=buildGeometry(opt);
[T]=linearConstr(unitCell,extrudedUnitCell,opt);
u0 = zeros(1, size(T,2)); % starting point
nsamples = 1e3;
nrepeats = 6;
func = 'breathing_gradient';

switch func
    case 'mh'
        [smpl, acceptance, energy, deformationV, maxTheta] = ...
            mh(opt, u0, nsamples);
        mh_plot(nsamples, energy, deformationV, opt, acceptance, maxTheta)
    case 'breathing'
        [smpl, acceptance, energy, deformationV, maxTheta] = ...
            breathing(opt, u0, nsamples, nrepeats);
        mh_plot(nsamples*nrepeats, energy, deformationV, opt, ...
            acceptance, maxTheta)
    case 'annealing'
        [smpl, acceptance, energy, deformationV, maxTheta] = ...
            annealing(opt, u0, nsamples);
        mh_plot(nsamples, energy, deformationV, opt, acceptance, maxTheta)
    case 'breathing_gradient'
        [smpl1, acceptance1, energy1, deformationV1, maxTheta1] = ...
            breathing(opt, u0, nsamples, nrepeats);
        mh_plot(nsamples*nrepeats, energy1, deformationV1, opt, ...
            acceptance1, maxTheta1)
        
        % start optimsations
        [~,idx] = max(maxTheta1);
        u0 = deformationV1(idx, :)';
        opt.Khinge = 1;
        
        % feed to fmincon
        opt.date = datestr(now, 'mmm-dd-yyyy');
        opt.time = datestr(now,'HH-MM-SS');
        [energy2, deformationV2] = grad_descent_fmincon(opt,u0);
        close all
        mh_plot(length(energy2), energy2, deformationV2, opt);
        
        % feed to naive gradient descent
        max_iter = 1e4;
        opt.date = datestr(now, 'mmm-dd-yyyy');
        opt.time = datestr(now,'HH-MM-SS');
        [energy2, deformationV2] = grad_descent(opt,u0,max_iter);
        close all;
        mh_plot(length(energy2), energy2, deformationV2, opt);
    otherwise
        error('Function not defined')
end