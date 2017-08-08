function [nodes, edges] = regraphNodes(nodes, edges)


% Re-graph 
nodePos = nodes;
nodeEdges = edges;

nNodes = size(nodePos,1);
nEdges = size(nodeEdges,1);
if ~exist('nodeDiam')
    nodeDiam = zeros(nNodes,1);
end

hxy = 2;
hz = 2;

nNodesUnique = 1;
nodeMap = zeros(nNodes,1);
nodeUnique = zeros(nNodes,1);

nodeMap(1) = 1;
nodePosNew = nodePos(1,:);
nodeUnique(1) = 1;

hwait = waitbar(0,'Regraphing nodes...');   %modify this if using parfor

for ii=2:nNodes
   if isequal(rem(ii,1000),0)
   waitbar(ii/nNodes,hwait);   %updateing waitbar takes a long time. 
   end
    pos = nodePos(ii,:);
    lst = find(pos(1)>=(nodePosNew(:,1)-hxy) & pos(1)<=(nodePosNew(:,1)+hxy) & ...
        pos(2)>=(nodePosNew(:,2)-hxy) & pos(2)<=(nodePosNew(:,2)+hxy) & ...
        pos(3)>=(nodePosNew(:,3)-hz) & pos(3)<=(nodePosNew(:,3)+hz) );
    if isempty(lst)
        nNodesUnique = nNodesUnique+1;

        nodeMap(ii) = nNodesUnique;
        nodeUnique(ii) = 1;
        
        nodePosNew(nNodesUnique,:) = pos;
        nodeDiamNew(nNodesUnique) = nodeDiam(ii);
    else
        if length(lst)>1
            clear d
            for iLst=1:length(lst)
                d(iLst) = norm(pos-nodePosNew(lst(iLst),:));
            end
            [foo, closestNode] = min(d);
        else
            closestNode = 1;
        end
        nodeMap(ii) = lst(closestNode);
    end
end
close(hwait);

nodeEdgesNew = nodeMap(nodeEdges);
nodeEdges = nodeEdgesNew;

%%%%%%%%%%%%%%
% prune edges - still need to handle small loops
% point edges
nodeEdges = nodeEdges(find(nodeEdges(:,1)~=nodeEdges(:,2)),:);
% redundant edges
sE = cell(size(nodeEdges,1),1);

for ii=1:length(nodeEdges)
    if nodeEdges(ii,1)<nodeEdges(ii,2)
        sE{ii} = sprintf('%05d%05d',nodeEdges(ii,1),nodeEdges(ii,2));
    else
        sE{ii} = sprintf('%05d%05d',nodeEdges(ii,2),nodeEdges(ii,1));
    end
end

[b,i,j]=unique(sE);
nodeEdges = nodeEdges(sort(i),:);

nodes = nodePosNew;
edges = nodeEdges;

disp(sprintf('Regraph reduced %d nodes to %d, and %d loops to %d',nNodes,size(nodes,1),...
    nEdges-nNodes+1,size(edges,1)-size(nodes,1)+1))

