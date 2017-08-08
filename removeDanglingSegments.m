function [seg, nodePos, nodeEdges] = removeDanglingSegments( seg, nodePos, minNvertices )


lstPrune = [];
for ii=1:length(seg)
    if seg(ii).nv < minNvertices 
        % check if dangling segment
        f1 = 0;
        f2 = 0;
        for jj=1:length(seg)
            if jj~=ii
                if ~isempty(find(seg(ii).v(1)==seg(jj).v))
                    f1 = 1;
                end
                if ~isempty(find(seg(ii).v(end)==seg(jj).v))
                    f2 = 1;
                end
            end
        end
        
        if f1==0 | f2==0
            lstPrune(end+1) = ii;
        end
    end
end
seg(lstPrune) = [];

% create nodeEdges
if ~isempty(seg)
    nodeEdges = seg(1).e;
    for ii=2:length(seg)
        nodeEdges = [nodeEdges; seg(ii).e];
    end
else
    nodeEdges = [];
end

% prune nodes and correct edges
nNodesPruned = 0;
lstNodesPruned = [];
for iN = 1:size(nodePos,1)
    lst = find(nodeEdges==iN);
    if length(lst)==0
        nNodesPruned = nNodesPruned + 1;
        lstNodesPruned(end+1) = iN;
    else
        nodeEdges(lst) = iN - nNodesPruned;
    end
end
nodePos(lstNodesPruned,:) = [];

% update seg structure
[im, seg] = findSegmentGroups( nodePos, nodeEdges );

