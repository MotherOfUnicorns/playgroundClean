function saveCSV(opt, header, fn_name, varargin)
% saveCSV(opt, header, fn_name, varargin)
%
% Dumps generated data into several CSV files.
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% opt
% header
% fn_name  - name of the function calling saveCSV.m
% varargin - must be in the form of [name1],[array1], [name2],[array2], ...
%            e.g., 'smpl', smpl, 'acceptance', acceptance
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% One csv file for each input variable
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on May 22, 2017
% yun

% create folder if it doesn't exist yet
date = opt.date;
time = opt.time;
nameFolder = [pwd, '/Results/', date, '/data/'];
if exist(nameFolder, 'dir')==0
    mkdir(nameFolder);
end

% check inputs
if nargin == 3 || mod(nargin, 2) == 0
    error('Oops, wrong inputs')
end

stuff = struct();
for ct = 1:(nargin-3)/2
    [stuff(:).(varargin{2*ct-1})] = deal(varargin{2*ct});
end

% save all files
fields = fieldnames(stuff);
for ct = 1:length(fields)
    file_name = strcat(opt.template, '_', time, '_', fn_name, '_',...
        fields(ct), '.csv');
    data = stuff.(fields{ct});
    full_name = char(fullfile(nameFolder, file_name));
    dlmwrite(full_name, header, 'delimiter', '', '-append')
    dlmwrite(full_name, data, '-append');
end