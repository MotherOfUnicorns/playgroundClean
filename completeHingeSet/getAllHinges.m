function getAllHinges(G, opt)
% get the complete list of hinge sets, and dump them in a csv file
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% G - a directed graph created from the options
% opt - options
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% A .csv file of the complete list of hinge sets, where each row is one
% such hinge set
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Jun 05, 2017
% yun


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