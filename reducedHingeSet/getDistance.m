function dis = getDistance(G, flavourDict)
% Returns the ``distance'' matrix with flavours
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% G - an digraph object
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% dis - a distance matrix
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Jun 07, 2017
% yun



% initialising...
H = height(G.Nodes);
dis = zeros(H);

% loop through everything, because we need to
for ii = 1:H
    for jj = 1:H
        % get all possible paths between two nodes *ii* and *jj*
        pathsiijj = pathbetweennodes(G.adjacency, ii, jj);
        % get the lengths of each of these paths
        pathsLeniijj = zeros(length(pathsiijj), 1);
        for pCt = 1:length(pathsiijj)
            pathsLeniijj(pCt) = length(pathsiijj{pCt});
        end
        
        % get all paths with the same minimum lengths, and only consider
        % these shortest ones for our new definition of ``distance''
        minPathLen = min(pathsLeniijj);
        minPathsIdx = find(minPathLen==pathsLeniijj);
        
        % get the ``distance'' from node *ii* to node *jj*
        % by looping through all minimum length paths
        pathDis = inf; % initialise with infinity
        for mCt = 1:length(minPathsIdx)
            currentPath = pathsiijj{minPathsIdx(mCt)};
            currentPathDisStr = '';
            
            % get the string of ``distance''
            for lCt = 1:length(currentPath)
            	currentPathDisStr = strcat(currentPathDisStr, ...
                    flavourDict(G.Nodes{currentPath(lCt),'Type'}{1}));
            end
            
            % convert to number and take the smallest value of the
            % ``distance'' and the flipped ``distance''
            currentPathDis = min(str2num(currentPathDisStr),...
                str2num(flip(currentPathDisStr)));
            if currentPathDis < pathDis
                pathDis = currentPathDis;
            end
        end

        % update distance matrix
        dis(ii,jj) = pathDis;
    end
end


% update distance matrix to make it symmetrical
for ii = 1:H
    for jj = ii:H
        dis(ii,jj) = min(dis(ii,jj), dis(jj,ii));
        dis(jj,ii) = dis(ii,jj);
    end
end
        