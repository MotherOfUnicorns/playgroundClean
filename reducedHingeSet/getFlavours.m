function [flavourDict, flavourTypes, flavourNum] = getFlavours(G)
% Returns the info about flavours of a graph G
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% G -  a directed graph
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% flavourTypes - a cell array containing the types of flavours
% flavourNum   - number of different flavours
% flavourDict  - a container map that maps different flavourTypes to
%                different numeric values. E.g., in truncated tetrahedron,
%                flavourDict('3-6') = 1, flavourDict('6-6') = 2
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Jun 07, 2017
% yun



flavourTypes = unique(G.Nodes.Type);
flavourNum = length(flavourTypes);

flavourDict = containers.Map;
for ii = 1:flavourNum
    flavourDict(flavourTypes{ii}) = num2str(ii);
end