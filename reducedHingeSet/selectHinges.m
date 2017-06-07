function [hinges] = ...
    selectHinges(G, dis, flavourTypes, flavourNum, N, hingeSetsPrev)
% select *N* number of hinges from graph *G* to actuate, the flavoured
% version
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% G - graph object
% dis - matrix containing the minimum distance between any two nodes
% N - number of hinges to be actuated
% hingeSetsPrev - the complete hinge sets where each set contains (N-1) hinges
% opt - options
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% hinges - an x*N array, where each row is a set of hinges to be actuated
%          at the same time
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Jun 07, 2017
% yun

% input validation
if N < 0 || N > size(G.Nodes,1)/2
    hinges = [];
    warning('You are selecting the wrong number of hinges to actuate!')

elseif N == 0
    hinges = [];
    
elseif N == 1
    % actuate different types of hinges
    hinges = zeros(flavourNum, 1);
    for ct = 1:flavourNum
        idx = find(1==strcmp(G.Nodes.Type, flavourTypes(ct)), 1);
        hinges(ct) = str2num(G.Nodes{idx,'Name'}{1});
    end
    
    
else% N >= 2
    % initialising...
    hinges = [];
    currentIdx = 1;
    cache.distance = []; % the matrices representing each selection
    cache.distanceEig = []; % their corresponding eigenvalues    

    
    % loop through each hinge set
    for ct = 1:size(hingeSetsPrev, 1)
        hingeSet = hingeSetsPrev(ct, :);
        hingeSetIdx = getHingeIdx(hingeSet, G);
        cols = dis(:, hingeSetIdx);
        
        % first, split the possible candidates into cases of different
        % flavours
        for flavourCt = 1:flavourNum
            candidateIdxs = find(1==strcmp(G.Nodes.Type, ...
                flavourTypes{flavourCt}));
            
            % then start comparing matrices
            for nodeCt = 1:length(candidateIdxs)
                newNodeIdx = candidateIdxs(nodeCt);
                if newNodeIdx <= hingeSetIdx(end)
                    % avoid double counting
                    continue;
                end
                
                newNode = str2num(G.Nodes{newNodeIdx, 'Name'}{1});
                newSet = [hingeSet, newNode];
                newSetIdx = [hingeSetIdx, newNodeIdx];
                newDisMat = dis(newSetIdx, newSetIdx);
                newDisEig = round(1000*sort(eig(newDisMat)')) / 1000; % round off 
                
                % if newDisEig is not in cache, add current solution
                if isempty(cache.distanceEig) || ...
                        ~ismember(newDisEig, cache.distanceEig, 'rows')
                    hinges(currentIdx, :) = newSet;
                    cache.distance(:,:,currentIdx) = newDisMat;
                    cache.distanceEig(currentIdx, :) = newDisEig;
                    currentIdx = currentIdx + 1;
                end        
                    
            end
            
        end
        
    end    

end



function hingeSetIdx = getHingeIdx(hingeSet, G)
% gets the index of hinges in graph, as numbered in *outputResults()*
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% hingeSet - an array of integers, which represent the names of hinges
% G -  a directed graph object
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% hingeSetIdx - an array of integers of the same size as hingeSet, which
%               contains the indeces of the hinges in the graph G
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Apr 10, 2017
% yun

hingeSetIdx = zeros(size(hingeSet));
for ct = 1:length(hingeSetIdx)
    hingeSetIdx(ct) = find(strcmp(G.Nodes.Properties.RowNames,...
                      num2str(hingeSet(ct))));
end