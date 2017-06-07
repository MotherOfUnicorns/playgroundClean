function getAllHinges(G, dis, flavourTypes, flavourNum, ...
    unitCell, extrudedUnitCell, opt, plotSelection)
% get the complete list of hinge sets, and dump them in a csv file
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% G - a directed graph created from the options
% flavourTypes
% flavourNum
% opt - options
% plotSelection - boolean for plotting selected hinges
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% A .csv file of the complete list of hinge sets, where each row is one
% such hinge set
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Jun 05, 2017
% yun

% input validation
% generate plots of selected hinges only when asked
if ~exist('plotSelection', 'var')
    plotSelection = false;
end

% create folder if it doesn't exist
folderName = strcat(pwd, '/hingeList_reduced/');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
fileName = strcat(folderName, opt.template, '.csv');
if exist(fileName, 'file')
    delete(fileName) % always start with new file
end

% name of all hinges
allHinges = zeros(1, height(G.Nodes));
for ct = 1:height(G.Nodes)
    allHinges(ct) = str2num(G.Nodes.Properties.RowNames{ct});
end

% starts generating all hignes
hingeSetsPrev = [];
for N = 1:(0.5 * height(G.Nodes))
    hingeSets = ...
        selectHinges(G, dis, flavourTypes, flavourNum, N, hingeSetsPrev);
    % and also initialise the complement selections
    cHingeSets = zeros(size(hingeSets, 1), height(G.Nodes)-N);
    
    for ct = 1:size(hingeSets, 1)
        % update the complement selections
        cHingeSets(ct,:) = setdiff(allHinges, hingeSets(ct,:));
        
        % write selections to file
        dlmwrite(fileName, hingeSets(ct,:), 'delimiter', ',', '-append')
        
        % plot and save figures when specified
        if plotSelection
            f = plotSelectedHinges(hingeSets(ct,:), unitCell, ...
                extrudedUnitCell, true);
            % set(f, 'visible', 'off');
            folderName = strcat(pwd, '/hingeList_reduced/', opt.template, '/', ...
                         num2str(N), '/');
            if ~exist(folderName, 'dir')
             mkdir(folderName)
            end
            saveas(f, strcat(folderName, mat2str(hingeSets(ct,:)), '.png'))
            close(f)
        end
    end
    
    
    for ct = 1:size(cHingeSets, 1)
        % write complement selections to file
        dlmwrite(fileName, cHingeSets(ct,:), 'delimiter', ',', '-append')
        
        % plot and save figures when specified
        if plotSelection
            f = plotSelectedHinges(cHingeSets(ct,:), unitCell, ...
                        extrudedUnitCell, true);
            % set(f, 'visible', 'off');
            folderName = strcat(pwd, '/hingeList_reduced/', opt.template, '/', ...
                        num2str(height(G.Nodes)-N), '/');
            if ~exist(folderName, 'dir')
                mkdir(folderName)
            end
            saveas(f, strcat(folderName,mat2str(cHingeSets(ct,:)), '.png'))
            close(f)
        end
    end

    % getting ready for the next loop
    hingeSetsPrev = hingeSets;
end