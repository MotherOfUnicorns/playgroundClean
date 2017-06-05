function plotConverged(opt)
% plotConverged(template)
% 
% plots the results that have converged to a different stable state other
% than the starting configuration. Also eliminates results that are badly
% scaled or ``not square enough''
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% opt - options
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% .png images of the selected results
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on May 22, 2017
% yun


bad = 0; % number of bad results
good = 0; % number of good results


[unitCell,extrudedUnitCell,opt]=buildGeometry(opt);
% keep a copy of extruded unit cell
extrudedUnitCell_original = extrudedUnitCell;


% get all files in directory that have exitFlag [1,1]
% i.e., those that converged correctly
fileFolder = strcat(pwd, '/Results/', opt.template, '/mat/');
finished = dir(fileFolder);
tol = .07;
for ct = 1:length(finished)
    if finished(ct).isdir
        disp('skip all directories...')
        continue;
    end
    
    % parse the file name to get back hinge set and exit flag
    fileName = finished(ct).name;
    parsedName = strsplit(fileName(1:end-4), '_');
    hingeSetStr = parsedName{1};
    exitFlagStr = parsedName{2};
    
    % convert to arrays of numbers
    hingeSetStr = strsplit(hingeSetStr(6:end), {' ', '[', ']'});
    exitFlagStr = strsplit(exitFlagStr(9:end-1), ' ');
    hingeSet = [];
    exitFlag = [0 0];
    for ii = 1:length(hingeSetStr)
        digit = hingeSetStr{ii};
        if ~strcmp('', digit)
            hingeSet = [hingeSet, str2num(digit)];
        end
    end
    exitFlag(1) = str2num(exitFlagStr{1});
    exitFlag(2) = str2num(exitFlagStr{2});
    
    % get rid of results that didn't converges
    if ~isequal(exitFlag, [1 1])
        disp('Discard result that does not converge...')
        bad = bad+1;
        continue;
    end
    
    % read data
    load(strcat(fileFolder, fileName));
    % update extrudedUnitCell
    extrudedUnitCell.node = extrudedUnitCell_original.node...
                          + result.deform(3).V;
    
    % get rid of results that went to NaN
    if isnan(extrudedUnitCell.node)
        disp('Discard result that does not converge correctly...')
        bad = bad + 1;
        continue;
    end
    
    % decide if it converged to a different stable state by checking if
    % there are hinges whose angles are close to pi or -pi
    u = result.deform(3).Ve;
    extrudedUnitCell.angleConstr = [];
    [~,~,~,~,~,theta] = Energy(u,extrudedUnitCell,opt);

    if sum(abs(theta)>pi)>0
        disp('Discard result that has faces crossing over...')
        bad = bad + 1;
        continue
    end
    if min(abs(abs(theta)-pi))>tol
        disp('Discard result that does not have zero angles...')
        bad=bad+1;
        continue;
    end
    
    % get rid of results that have not-very-square faces
    nextRound = false;
    for fCt = 1:length(extrudedUnitCell.face)
        if nextRound
            break;
        end
        face = extrudedUnitCell.face{fCt};
        diag1 = face([1 3]);
        diag2 = face([2 4]);
        diag1n1 = extrudedUnitCell.node(diag1(1),:);
        diag1n2 = extrudedUnitCell.node(diag1(2),:);
        diag2n1 = extrudedUnitCell.node(diag2(1),:);
        diag2n2 = extrudedUnitCell.node(diag2(2),:);
        diag1len = norm(diag1n1-diag1n2, 2);
        diag2len = norm(diag2n1-diag2n2, 2);
        if abs(diag1len-diag2len) > tol
            nextRound = true;
        end
    end
    if nextRound
        disp('Discard result that have faces not square enough...')
        bad = bad+1
        continue;
    end
    
    
    % plot stable state and save
    newResult.numMode = 1;
    newResult.deform.V = result.deform(3).V;
    newResult.deform.Ve = result.deform(3).Ve;
    outputResults(unitCell,extrudedUnitCell,newResult,opt);
    fig = gcf;
    set(fig, 'visible', 'off');

    imgFolder = strcat(pwd, '/Results/', opt.template, '/img/');
    if ~exist(imgFolder, 'dir')
        mkdir(imgFolder)
    end
    saveas(fig, strcat(imgFolder, fileName(1:end-3), 'png'))
    close all;
    
    good = good + 1
end