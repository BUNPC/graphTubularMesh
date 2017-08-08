function [im, segs] = findSegmentGroups( nodePos, nodeEdges )

%%%%%%%%%%
% TO DO:
% correct eLen (segLength) to be in um (non-cubical voxels problem) -  easy
% correct segDiam to be in um (non-cubical voxels problem) - more
% complicated - we already have calculated diameters in imView in voxels
%%%%%%%%%%


hwait = waitbar(0,'Getting Segment and Group Info');

nB = numNodeEdges( nodePos, nodeEdges );

nN = size(nodePos,1);

hxy = 1;%im.Hvox(1);
hz  = 1;%im.Hvox(3);

% create nodePos_um variable, to have positions in um
nodePos_um = nodePos; %this is a simple fix but we also need to
                   %update im.nodePos_um whenever we update im.nodePos
nodePos_um(:,1:2) = nodePos(:,1:2).*hxy;
nodePos_um(:,3) = nodePos(:,3).*hz;

lst3p = find(nB>=3 | nB==1);

edgeSegN = zeros(size(nodeEdges,1),1);
nodeSegN = zeros(size(nodePos,1),1);
nodeGrp = zeros(nN,1);
nSeg = 1;
nGrp = 0;
segLen = [];
segLen_um = [];
segNedges = [];
segEndNodes = [];
for ii = 1:length(lst3p)
    if isequal(rem(ii,100),0)
        waitbar(ii/length(lst3p),hwait);
    end
    iN = lst3p(ii);
    if nodeGrp(iN)==0
        nGrp = nGrp + 1;
        nodeGrp(iN) = nGrp;
    end
    nGrpN = nodeGrp(iN); %nGrp;
    lstE = find(nodeEdges(:,1)==iN | nodeEdges(:,2)==iN);
    iNstart = iN;
    for jj = 1:length(lstE)
        iN = lst3p(ii);
        eIdx =lstE(jj);
        % correct this to be in um (non-cubical voxels)
        eLen = sum(diff( nodePos(squeeze(nodeEdges(eIdx,:)),:), 1, 1).^2).^0.5;
        eLen_um = sum(diff( nodePos_um(squeeze(nodeEdges(eIdx,:)),:), 1, 1).^2).^0.5;      
        iN = setdiff(unique( nodeEdges(eIdx,:)), iN);
        if nodeGrp(iN)==0
            nodeGrp(iN) = nGrpN;
        elseif nodeGrp(iN) < nGrpN
            lst = find(nodeGrp==nGrpN);
            nodeGrp(lst) = nodeGrp(iN);
            nGrpN = nodeGrp(iN);
        elseif nodeGrp(iN) > nGrpN
            % I think this is fine
            lst = find(nodeGrp==nodeGrp(iN));
            nodeGrp(lst) = nGrpN;
%            if nB(iN)<3
%                error('We should never get here');
%            end
        end
        
%         if nodeSegN(iN)==0  % remove 6/3/09 since it appears below
%             nodeSegN(iN) = nSeg; % and should resolve an issue in
%                                  % nodeGrps()
%             nodeSegN(iNstart) = nSeg;            
%         end
        if nodeSegN(iN)==0 % added back 8/3/09 because if had nodeSegN=0 one node from bifurcation
            nodeSegN(iN) = nSeg;
        end

        nE = 1;
        nLst = [];
        if edgeSegN(eIdx)==0
            while nB(iN)==2
                nLst(end+1) = iN;
                edgeSegN(eIdx) = nSeg;
                lstE2 = find(nodeEdges(:,1)==iN | nodeEdges(:,2)==iN);
                eIdx = setdiff(lstE2,eIdx);
                % correct this to be in um (non-cubical voxels)
                eLen = eLen + sum(diff( nodePos(squeeze(nodeEdges(eIdx,:)),:), 1, 1).^2).^0.5;
                eLen_um = eLen_um + sum(diff( nodePos_um(squeeze(nodeEdges(eIdx,:)),:), 1, 1).^2).^0.5;
                iN = setdiff(unique(nodeEdges(eIdx,:)), iN);
                if nodeGrp(iN)==0
                    nodeGrp(iN) = nGrpN;
                elseif nodeGrp(iN) < nGrpN 
                    lst = find(nodeGrp==nGrpN);
                    nodeGrp(lst) = nodeGrp(iN);
                    nGrpN = nodeGrp(iN);
                elseif nodeGrp(iN) > nGrpN
                    lst = find(nodeGrp==nodeGrp(iN));
                    nodeGrp(lst) = nGrpN;
                end
                if nodeSegN(iN)==0
                    nodeSegN(iN) = nSeg;
                end
                nE = nE + 1;
            end
            nodeSegN(iNstart) = nSeg; % added 6/3/09 when deleted above
            nodeSegN(iN) = nSeg; % added 6/3/09 to resolve issue with nodeGrps
                                 % i hope this doesn't cause trouble
            edgeSegN(eIdx) = nSeg;
            segNedges(nSeg) = nE;
            segLen(nSeg) = eLen;
            segLen_um(nSeg) = eLen_um;
%            if ~isempty(nLst)
%                segDiam(nSeg) = median(im.nodeDiam(nLst));
%                segVesType(nSeg) = median(im.nodeType(nLst));
%            else
%                segDiam(nSeg) = mean(im.nodeDiam(im.nodeEdges(eIdx,:)));
%                segVesType(nSeg) = max(im.nodeType(im.nodeEdges(eIdx,:)));
%            end
            kk = find(lst3p==iN);
            if ~isempty(kk) 
                segEndNodes(end+1,:) = lst3p([ii kk]);
            else
                error('We should not get here')
            end
            nSeg = nSeg + 1;
        end
    end
end
close(hwait)

% remove groups with zero nodes
nGrpN = 0;
for ii=1:nGrp
    lst = find(nodeGrp==ii);
    if length(lst)>0
        nGrpN = nGrpN + 1;
        nodeGrp(lst) = nGrpN;
    end
end
nGrp = nGrpN;


% save stats
im.nodeGrp = nodeGrp;
im.nodeSegN = nodeSegN;
im.segNedges = segNedges(:);
im.segLen = segLen(:);
im.segLen_um = segLen_um(:);
%im.segDiam = segDiam;
%im.segVesType = segVesType;
im.edgeSegN = edgeSegN;
im.segEndNodes = segEndNodes;
%im.segNodeMap = lst3p;
%im.nB = nB;
im.segPos = squeeze(mean(reshape(nodePos(segEndNodes,:),[2 length(segLen) 3]),1));

grpNnodes = [];
disp( sprintf('\n\n# Groups = %d\nGrp\t\tE1\t\tE2\t\tE3\t\tE4+\t\tTotal',nGrp) )
for ii=1:nGrp
    lst = find(nodeGrp==ii);
    grpNnodes(ii) = length(lst);
    im.Stats.grpE1234(ii,1) = length(find(nB(lst)==1));
    im.Stats.grpE1234(ii,2) = length(find(nB(lst)==2));
    im.Stats.grpE1234(ii,3) = length(find(nB(lst)==3));
    im.Stats.grpE1234(ii,4) = length(find(nB(lst)>3));
    im.Stats.grpE1234(ii,5) = sum(im.Stats.grpE1234(ii,1:4));
    disp(sprintf('%d\t\t',[ii im.Stats.grpE1234(ii,:)]))
end
im.Stats.grpNnodes = grpNnodes;
 

segs = [];
for iSeg = 1:length(segLen)
    lstE = find(edgeSegN==iSeg);
    e = nodeEdges(lstE,:);
    segs(iSeg).e = e;
    
    % needs to be in order of connected nodes from one end of segment to the other
    iN = segEndNodes(iSeg,1);
    lstV = iN;
    lst1 = find(e(:,1)==iN | e(:,2)==iN);
    iN2 = setdiff( e(lst1,:), lstV );
    while ~isempty(iN2)
        lstV(end+1) = iN2;
        lst1 = find(e(:,1)==iN2 | e(:,2)==iN2);
        iN2 = setdiff( e(lst1,:), lstV );
    end
    segs(iSeg).v = lstV; 
    
    segs(iSeg).ne = size( segs(iSeg).e, 1 );
    segs(iSeg).nv = length( segs(iSeg).v );
end

