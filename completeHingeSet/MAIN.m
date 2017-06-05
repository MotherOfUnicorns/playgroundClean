% generates the list of all hinge sets and saves them in a .csv file,
% reads from then file and do 2-step optimisation, and finally
% plots the converged meaningful results
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Jun 05, 2017
% yun

template = 'tetrahedron';
opt=initOpt('inputType','individual',...
            'template',template,...
            'plot','result',...
            'interval',1,'saveFig','off','periodic','on',...
            'figDPI',300,...
            'saveMovie', 'on', 'safeMovieAntiAlias', 0,...
            'constrFace','off','constrEdge','off',...
            'Khinge',0.0005,'Kedge',1,'Kface',1,'KtargetAngle',1,...
            'constAnglePerc',0.99);
[unitCell, extrudedUnitCell, opt] = buildGeometry(opt);
[G, opt] = buildGraph(unitCell, extrudedUnitCell, opt);

%% first generate the list of all hinge sets
getAllHinges(G, opt);


%% then perform the 2-step optimisation and save the resutls in .mat files
fileName = strcat(pwd, '/hingeList/', template, '.csv');
% get the list of hinges need to actuate
hingeList = dlmread(fileName);

for ct = 1:size(hingeList, 1)
    row = hingeList(ct, :);
    close all
    hinges = row(0~=row); % trim off the trailing zeros
    
    [result, extrudedUnitCell, opt] = ...
        actuateHinges(hinges, unitCell, extrudedUnitCell, opt);
end


%% finally, plot the meaningful results
plotConverged(opt)
