function [result, extrudedUnitCell, opt] = actuateHinges(hingeList, opt)
% [result, extrudedUnitCell, opt] = actuateHinges(hingeList, opt)
% 
% Two-step optimiation with hinges specified in the hingeList closed.
% Also saves the the results and exitFlags in a .mat file, (and after
% uncommenting the final section, outputs an image of the final state).
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% hingeList - a list of hinges to be actuated at the same time
% opt       - options
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% result           - final results from convergence
% extrudedUnitCell - updated extruded unit cell
% opt              - updated options
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


% create folder if it doesn't exist
folderName = strcat(pwd, '/Results/', opt.template, '/mat/');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
fileName = strcat(folderName, 'hinge', mat2str(hingeList), '_exitflg',...
    mat2str([opt.exitFlag1, opt.exitFlag2]), '.mat');
save(fileName, 'result');

% % uncomment to output the final image
%     if opt.exitFlag1 == 1 && opt.exitFlag2 == 1
%         resultF.numMode = 1;
%         resultF.deform.V = result.deform(3).V;
%         resultF.deform.Ve = result.deform(3).Ve;
%         outputResults(opt.unitCell, extrudedUnitCell, resultF, opt)
%     end