% a script to save .fig as .png in directory

if ~exist('newnames\', 'dir')
    mkdir(pwd, 'newnames\')
end

% read all .fig files in current directory
files = dir(pwd);
names = {files.name};
newnames = {'mh', 'fmincon', 'grad'};

% first loop for figs
figname_ct = 0;
nr_ct = 1;
for ct = 1:length(names)
    name = names{ct};
    
    if strfind(name, '.fig')
        f = openfig(name);
        figname_ct = mod(figname_ct, 3);
        figname = [pwd '\newnames\' num2str(nr_ct) '_' ...
            newnames{figname_ct+1} '.png'];
        saveas(f, figname, 'png');
        
        figname_ct = figname_ct + 1;
        if mod(figname_ct, 3) == 0
            nr_ct = nr_ct + 1;
        end
    end

    close all
end



% second loop for gifs
gifname_ct = 0;
nr_ct = 1;
for ct = 1:length(names)
    name = names{ct};
    
    if strfind(name, '.gif')
        
        gifname_ct = mod(gifname_ct, 3);
        gifname = [pwd '\newnames\' num2str(nr_ct) '_' ...
            newnames{gifname_ct+1} '.gif'];
        copyfile(name, gifname)
        
        gifname_ct = gifname_ct + 1;
        if  mod(gifname_ct, 3) == 0
            nr_ct = nr_ct + 1;
        end
    end

    close all
end