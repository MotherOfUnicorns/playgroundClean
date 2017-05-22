function [result, extrudedUnitCell, opt] = actuateHinges(hingeList, opt)
% close hinges specified in the list
% also saves the hinges actuated the resulting exitflag, and the results in
% a .mat file, and output an image of the final stage
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% hingeList - a list of hinges to be actuated at the same time
% opt       - options
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Apr 04, 2017
% yun


if isempty(hingeList)
    result = [];
    extrudedUnitCell = opt.extrudedUnitCell;
    return
end

newList = [hingeList(:), -pi*0.985 * ones(length(hingeList), 1)];
opt.angleConstrFinal(1).val = newList;
opt.options=optimoptions('fmincon','GradConstr','on','GradObj','on',...
    'tolfun',1e-5','tolx',1e-9,'tolcon',1e-5,'Display','off',...
    'DerivativeCheck','off','maxfunevals',100000);

[unitCell,extrudedUnitCell,opt]=buildGeometry(opt); % Bas' magic line
opt.unitCell = unitCell;
opt.extrudedUnitCell = extrudedUnitCell;
[result, extrudedUnitCell, opt] = ...
    findDeformation(unitCell,extrudedUnitCell,opt);


date = datestr(now, 'mmm-dd-yyyy');
% create folder if it doesn't exist
folderName = strcat(pwd, '/Results/', opt.template, '/', date, '/');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
fileName = strcat(folderName, 'hinge', mat2str(hingeList), '_exitflg',...
    mat2str([opt.exitFlag1, opt.exitFlag2]), '.mat');
save(fileName, 'result');

% % output the final image
%     if opt.exitFlag1 == 1 && opt.exitFlag2 == 1
%         resultF.numMode = 1;
%         resultF.deform.V = result.deform(3).V;
%         resultF.deform.Ve = result.deform(3).Ve;
%         outputResults(opt.unitCell, extrudedUnitCell, resultF, opt)
%     end