function f = plotSelectedHinges(hingeList, unitCell, extrudedUnitCell, ...
             orderedNames)
% function plotSelectedHinges(geometry, hingeList)
% highlights the set of selected hinges in a given polyhedron, for better
% visualisation
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% hingelist - a list containing the index of the selected hinges
% unitCell
% extrudedUnitCell
% orderedNames - boolean indicating whether the hinge names are printed as
%                the original order as in 
%                extrudedUnitCell.nodeHingeEx --> orderedNames = false
%                or in ordered version --> orderedNames = true
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% f - a graph object
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Jun 07, 2017
% yun


if ~exist('orderedNames', 'var')
    orderedNames = true;
end

vertices = unitCell.Polyhedron.node;

% plot transparent polyhedra with patch()
faces = unitCell.Polyhedron.face;
figure; hold on;
xlim([-1 1]); ylim([-1 1]); zlim([-1 1]); view([20 10]); axis('equal')
axis('off')
xlabel('x'); ylabel('y'); zlabel('z')
if 1==length(hingeList)
    hingeStr = num2str(hingeList);
else
    hingeStr = mat2str(hingeList);
end

title(strcat('Hinges selected: ',hingeStr(2:(end-1))))
for ct = 1:length(faces)
    patch('Vertices',vertices, 'Faces',faces{ct},...
          'FaceColor',[.73 .86 .95], 'FaceAlpha',0.3, ...
          'EdgeColor',[.2 .2 .2], 'lineWidth',0.5);
end
light('Position',[1 -1 1],'Style','local')


% find corresponding nodes on the hingeList and highlight those edges
for ed = 1:length(hingeList)
    endNodes = extrudedUnitCell.nodeHingeEx(hingeList(ed), 1:2);
    X = vertices(endNodes, 1);
    Y = vertices(endNodes, 2);
    Z = vertices(endNodes, 3);
    line(X, Y, Z, 'lineStyle',':', 'color', [.01 .23 .45], 'LineWidth', 8)
end


% add lables for the edges
for ct = 1:size(unitCell.Polyhedron.edge, 1)
    endNodes = unitCell.Polyhedron.edge(ct, :);
    
    if orderedNames
        hingeName = num2str(ct);
    else
        [~, hinge] = ismember(endNodes,...
                     extrudedUnitCell.nodeHingeEx(:,1:2),'rows');
        hingeName = num2str(hinge);
    end
    
    textPos = 0.5 * (vertices(endNodes(1), :) + vertices(endNodes(2), :));
    text(textPos(1)+.01, textPos(2)+.01, textPos(3)-.04, hingeName, ...
        'FontSize', 26, 'FontWeight', 'bold', 'color', [.2 .2 .2]);
end

f = gcf;