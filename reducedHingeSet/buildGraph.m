function [G, flavourTypes] = buildGraph(unitCell, extrudedUnitCell)
% build directed graph from polyhedra
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% unitCell
% extrudedUnitCell
% opt - options
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% G   - a (directed) graph object
% opt - updated options
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% EXAMPLE
% opt=initOpt('inputType','individual',...
%             'template','truncated tetrahedron',...
%             'plot','result',...
%             'interval',1,'saveFig','off','periodic','off',...
%             'constrFace','off','constrEdge','off',...
%             'Khinge',0.0005,'Kedge',1,'Kface',1,'KtargetAngle',1,...
%             'constAnglePerc',0.99);
% [unitCell, extrudedUnitCell, opt] = buildGeometry(opt);
% [G, opt] = buildGraph(unitCell, extrudedUnitCell, opt);
% plot(G)
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Jun 05, 2017 
% yun



% initialising...
numEdges = size(unitCell.Polyhedron.edge, 1);
G_adjacency = zeros(numEdges); % initialise adjacency matrix


% ========== uncomment to plot info ==========
% outputResults(unitCell, extrudedUnitCell, [], opt);


% center of the polyhedron, used for determining directions in graph *G*
center = mean(unitCell.Polyhedron.node, 1);


% determine how two hinges are connected, and generate adjacency matrix
hingeNames = cell(numEdges, 1); % corresponds to names in extrudedUnitCell
hingeTypes = cell(numEdges, 1);
isClosed = cell(numEdges, 1);
hingeNodes = cell(numEdges, 1);
for ii = 1:numEdges
    n1 = unitCell.Polyhedron.edge(ii,1);
    n2 = unitCell.Polyhedron.edge(ii,2); % endnodes of the edge
    hingeNodes{ii} = sort([n1 n2]);
    % getting the edge numbering from nodeHingeEx (consistent with how
    % they're numbered in outputResults
    [~,idx1]=ismember([n1 n2],extrudedUnitCell.nodeHingeEx(:,1:2),'rows');
    [~,idx2]=ismember([n2 n1],extrudedUnitCell.nodeHingeEx(:,1:2),'rows');
    hingeNames{ii} = num2str(idx1+idx2);
    
    % getting the types of faces associated with current edge
    faceTypes = [0 0];
    for ff = 1:length(unitCell.Polyhedron.face)
        n_face = unitCell.Polyhedron.face{ff};
        if sum(find(n1==n_face)) && sum(find(n2==n_face))
            faceTypes(find(0==faceTypes, 1)) = ...
                length(unitCell.Polyhedron.face{ff});
        end
    end
	hingeTypes{ii} = [num2str(min(faceTypes)),'-',num2str(max(faceTypes))];
    isClosed{ii} = false;
    
    % here starts stuff for adjacency matrix
    for jj = 1:numEdges
        if ii == jj
            G_adjacency(ii,jj) = 0;
        else
            n_edge2 = unitCell.Polyhedron.edge(jj, :);
            
            % see if these two hinges share a common node
            comm_node = [n_edge2(n1==n_edge2), n_edge2(n2==n_edge2)];
            if ~isempty(comm_node)
                % determine the direction between these two hinges
                edgeV1 = - getHingeVector([n1,n2], comm_node, unitCell);
                edgeV2 = getHingeVector(n_edge2, comm_node, unitCell);
                normV = unitCell.Polyhedron.node(comm_node,:) - center;
                direction = dot(normV, cross(edgeV1,edgeV2));
                if direction < 0
                    G_adjacency(ii,jj) = 1;
                elseif direction > 0
                    G_adjacency(jj,ii) = 1;
                end
            end
        end
    end
end
G_adjacency = sparse(G_adjacency);

% finally, make a graph of the hinges
G = digraph(G_adjacency);

% assign names and other properties to the nodes of the graph (edges on the polyherdon)
G.Nodes.Name = hingeNames;
G.Nodes.Properties.RowNames = hingeNames;
G.Nodes.IsClosed = isClosed;
G.Nodes.Type = hingeTypes;
G.Nodes.HingeNodes = hingeNodes;

% opt.graph = G;
end


function [vec] = getHingeVector(nodes, comm_node, unitCell)
% returns the hinge in the form of a 3D vector, and always use *comm_node*
% as the starting point of the vector
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% nodes     - a 2*1 array, containing the two nodes of a hinge
% comm_node - the node in common between this hinge and the other one 
%             under consideration
% unitCell  - a unit cell
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% vec - a vector of this hinge, pointing from *comm_node* to the other node
%       in *nodes*
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Mar 20, 2017


comm_idx = find(comm_node==nodes);
if comm_idx == 1
    % do nothing
elseif comm_idx == 2
    nodes = flip(nodes);
else
    error('Oops, you''ve got the wrong nodes!')
end

coord1 = unitCell.Polyhedron.node(nodes(1), :);
coord2 = unitCell.Polyhedron.node(nodes(2), :);
vec = coord2 - coord1;
end