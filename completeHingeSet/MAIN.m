% reads hingeList from file, and does the actuation and stuff
% last modified on Mar 31, 2017
% yun



template = 'truncated tetrahedron';
fileName = strcat(pwd, '/hingeList/', template, '.csv');
hingeList = dlmread(fileName);


% for ct = 100:100
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

