function nodes = straightenNodes( nodes, edges, Mask, Ithresh )

% Move towards mean of neighboring nodes, i.e. straigthening
nN = size(nodes,1);
nodePos = nodes;
nLst = 1:size(nodePos,1);
hwait = waitbar(0,'Moving nodes towards mean of neighboring nodes');
for iii= 1:length(nLst) %1:nN
    ii = nLst(iii);
    waitbar(iii/length(nLst),hwait)
    eLst = find(edges(:,1)==ii | edges(:,2)==ii);
    nLst2 = setdiff(unique(edges(eLst,:)), ii);
    
    proceedFlag=1; % filter nodes on diamter thresshold
    proceedFlag2=1; % select nodes inside Z range
    proceedFlag3=1; % select nodes inside XY range
    
    if proceedFlag && proceedFlag2 && proceedFlag3
        if length(nLst2)>1
            pos0 = max(nodes(ii,:),1);
            posC = mean(nodes(nLst2,:),1);
            posN = pos0 + (posC-pos0) / max(norm(posC-pos0),1);
            %            if im.I(round(pos0(2)),round(pos0(1)),round(pos0(3)))>=Ithresh
            if Mask(round(posN(2)),round(posN(1)),round(posN(3)))>=Ithresh
                nodePos(ii,:) = posN;
            end
            %            else
            %                nodePos(ii,:) = posN;
            %            end
        end
    end
    
end
close(hwait)
nodes = nodePos;

