% generates the list of all hinge sets and saves them in a .csv file,
% reads from then file and do 2-step optimisation, and finally
% plots the converged meaningful results
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on May 22, 2017
% yun

template = 'tetrahedron';

%% first generate the list of all hinge sets
getAllHinges(template);


%% then perform the 2-step optimisation and save the resutls in .mat files
fileName = strcat(pwd, '/hingeList/', template, '.csv');
hingeList = dlmread(fileName);

for ct = 1:size(hingeList, 1)
    row = hingeList(ct, :);
    close all
    hinges = row(0~=row);%hingeList(ct,0~=hingeList(ct, :));
    
    opt=initOpt('inputType','individual',...
                'template',template,...
                'plot','result',...
                'saveFig','off','periodic','on',...
                'interval', 1, 'figDPI',300,...
                'saveMovie', 'on', 'safeMovieAntiAlias', 0,...
                'constrFace','off','constrEdge','off',...
                'Khinge',0.0005,'Kedge',1,'Kface',1,'KtargetAngle',0.5,...
                'constAnglePerc',0.99);
    [opt.unitCell, opt.extrudedUnitCell,~] = buildGeometry(opt);
    [result, extrudedUnitCell, opt] = actuateHinges(hinges, opt);
end


%% finally, plot the meaningful results
plotConverged(template, opt)
