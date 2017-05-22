function [hinges] = selectHinges(G, N)
% select *N* number of hinges from graph *G* to actuate -- a complete set,
% without considering symmetry
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% G - graph object
% N - number of hinges to be actuated
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% hinges - an x*N array, where each row is a set of hinges to be actuated
%          at the same time
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Mar 28, 2017
% yun


% input validation
if N < 0 || N > size(G.Nodes,1)
    hinges = [];
    warning('You are selecting the wrong number of hinges to actuate!')
else
    % the name of all hinges (in extrudedUnitCell)
    numHinges = size(G.Nodes, 1);
    allHinges = zeros(1, numHinges);
    for ct = 1:numHinges
        allHinges(ct) = str2num(G.Nodes.Name{ct});
    end
    
    hinges = combnk(allHinges, N); % all combinations
end