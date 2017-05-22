% leftoverHinges
% picks out the hinges that were left over during a previously failed
% simulation, and restart those again :'(
% last modified on Apr 18, 2017
% yun


template = 'truncated tetrahedron';
newFileName = strcat(pwd, '\hingeList\', template, '_unfinished.csv');


% read hinges in hingeList
hingeFile = strcat(pwd, '\hingeList\', template, '.csv');
allHinges = csvread(hingeFile);

% read file names of the finished ones
fileFolder = strcat(pwd, '\Results\', template, '\dataComplete');
finished = dir(fileFolder);
for ct = 1:length(finished)
    if finished(ct).isdir
        disp('skip all directories...')
        continue;
    end
    hingesCellarray = strsplit(finished(ct).name, '_');
    hingesCellarray = hingesCellarray{1}(7:(end-1));
    hingesCellarray = strsplit(hingesCellarray, ' ');
    
    hinges = zeros(1, size(allHinges,2));
    for ii = 1:size(hingesCellarray, 2)
        temp = str2num(hingesCellarray{ii});
        hinges(1,ii) = temp;
    end
    
    
    % delete hinges that are already worked on
    [flg, idx] = ismember(hinges, allHinges, 'rows');
    if flg
        allHinges = allHinges([1:(idx-1) (idx+1):end], :);
    end
end

% save the hinge sets that were not worked on
for ct = 1:size(allHinges, 1)
    hingeSet = allHinges(ct, 1);
    dlmwrite(newFileName, hingeSet(0~=hingeSet), ...
        'delimiter', ',', '-append');
end

% then work on those with MAIN again...