function getAllHinges(geometry)
% get the complete set of hinges, and dump them in a csv file
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% geometry - a string of the name of the geometry
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Mar 31, 2017
% yun


opt=initOpt('inputType','individual',...
            'template',geometry,...
            'plot','info',...
            'interval',1,'saveFig','off','periodic','on',...
            'figDPI',300,...
            'saveMovie', 'on', 'safeMovieAntiAlias', 0,...
            'constrFace','off','constrEdge','off',...
            'Khinge',0.0005,'Kedge',1,'Kface',1,'KtargetAngle',1,...
            'constAnglePerc',0.99);
[G, opt] = buildGraph(opt);

% create folder if it doesn't exist
folderName = strcat(pwd, '/hingeList/');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
fileName = strcat(folderName, opt.template, '.csv');
if exist(fileName, 'file')
    delete(fileName) % always start with new file
end

% starts generating all hignes
for N = 1:size(G.Nodes, 1)
    hinges = selectHinges(G, N);
    for ct = 1:size(hinges, 1)
        dlmwrite(fileName, hinges(ct,:), 'delimiter', ',', '-append')
    end
end