%%
% TO DO
% change nodeEdges to edges throughout

function graphTubularMesh( f, v, Mask, offset )

%%
% load mesh.mat output from vesSegment
% This has Mask, f, v

% load mesh.mat

hwait = waitbar(0,'Graphing...');

vLeftInd = ones(size(v,1),1);

nGrps = 0;
segAll = [];

while length(find(vLeftInd==1))>0
    
    %%
    % select a random vertex
    lst = find(vLeftInd==1);
    vSel = lst( ceil(rand(1)*length(lst)) );
    
    waitbar(1-length(lst)/length(vLeftInd), hwait, sprintf('Graphing... finding vertices for group %d',nGrps+1) );

    %%
    % find all vertices and edges connected to the selected vertex
    vTmp = vSel;
    vCon = vSel;
    fCon = [];
    
    c = 0;
    while ~isempty( vTmp )
        c = c + 1;
        vv = vTmp(end);
        vTmp(end) = [];
        [ir,ic] = find( f==vv );
        
        lstv = unique(f(ir,:));
        lstv2 = setdiff( lstv, vCon );
        vCon = [vCon; lstv2(:)];
        vTmp = [vTmp; lstv2(:)];
        
        ir2 = setdiff( ir, fCon);
        fCon = [fCon; ir2];
        
        if mod(c,100)==0
            [length(vCon) length(vTmp) length(fCon)]
        end
    end
    
    vLeftInd(vCon) = 0;
    
    %%
    % Calculate number of faces for each edge
    waitbar(1-length(lst)/length(vLeftInd), hwait, sprintf('Graphing... faces for each edge for group %d',nGrps+1) );

    Medges = spalloc(size(v,1),size(v,1),size(v,1)*10);
    
    for ii = 1:length(fCon)
        lstv = f(fCon(ii),:);
        Medges(lstv(1),lstv(2)) = Medges(lstv(1),lstv(2)) + 1;
        Medges(lstv(1),lstv(3)) = Medges(lstv(1),lstv(3)) + 1;
        Medges(lstv(2),lstv(3)) = Medges(lstv(2),lstv(3)) + 1;
        Medges(lstv(2),lstv(1)) = Medges(lstv(2),lstv(1)) + 1;
        Medges(lstv(3),lstv(1)) = Medges(lstv(3),lstv(1)) + 1;
        Medges(lstv(3),lstv(2)) = Medges(lstv(3),lstv(2)) + 1;
    end
    
    %%
    % find vertices of edges with only one face. And find the first connected group.
    waitbar(1-length(lst)/length(vLeftInd), hwait, sprintf('Graphing... group %d',nGrps+1) );

    [lstV1, ic] = find(Medges==1);
    lstV1 = unique( lstV1 );
    
    % find connected groups
    [nCGhold, lstVcgHold] = findConnectedGroups( lstV1, Medges );
    for ii=1:nCGhold
        lstVcgHold(ii).nvg_con = 0;
    end
    
    %%
    % select the first moveable group of vertices
    if length(lstVcgHold)>0
            nGrps = nGrps + 1;
            
            lstVmoveable = lstVcgHold(1).v;
            lstVcgHold(1) = [];
            
            % initialize some variables
            nvg = 0; % start with 0 graph nodes (number of vertices in the graph, nvg)
            gV = []; % positions of the graph vertices
            vConLeft = vCon; % connected mesh vertices left to be graphed
            fConLeft = fCon; % connected faces to be graphed
            camera.CP = []; % camera information
            
            nSeg = 1; % number of graphed segments
            seg = []; % graph segment structure
            seg(1).v = [];
            seg(1).e = []; % edges reference vertex indices in gV, not seg.v, since the edges linearly connect vertices in seg.v
            seg(1).nv = 0;
            seg(1).ne = 0;
            
            %%
            % Graph
            graphTubularMesh_sub
            
            [seg2, gV2, gE2] = removeDanglingSegments( seg, gV, 4 ); % this results in abandon nodes in gV
%            displayMeshTube( f,v, fCon, [], [], gV, seg2, camera, 0.1, 3 )
            
            [gV3, gE3] = fillNodes( gV2, gE2 );
            nB = numNodeEdges( gV3, gE3 );
            [im, seg3] = findSegmentGroups( gV3, gE3 );
%            displayMeshTube( f,v, fCon, [], [], gV3, seg3, camera, 0.1, 4 )
            
            segAll{nGrps}.seg = seg3;
            segAll{nGrps}.nodes = gV3;

    end
end

%%
% put all segs together
nSegs = 0;
nodes = [];
edges = [];
seg = [];

for iSegAll = 1:length(segAll)
    nNodes = size(nodes,1);
    nNodesNew = size(segAll{iSegAll}.nodes,1);
    if nNodesNew>0 & ~isempty(segAll{iSegAll}.seg)
        nodes(nNodes+[1:nNodesNew],1:3) = segAll{iSegAll}.nodes;
        for iSeg = 1:length(segAll{iSegAll}.seg)
            nSegs = nSegs + 1;
            seg(nSegs).e = segAll{iSegAll}.seg(iSeg).e + nNodes;
            seg(nSegs).v = segAll{iSegAll}.seg(iSeg).v + nNodes;
            
            nEdges = size(edges,1);
            nEdgesNew = size(seg(nSegs).e,1);
            edges(nEdges+[1:nEdgesNew],1:2) = seg(nSegs).e;
        end
    end
end

%displayMeshTube( f,v, 1:size(f,1), [], [], nodes, seg, camera, 0.1, 4 )

nodes0 = nodes;

%%
nodes = nodes0 + ones(size(nodes,1),1) * offset;

save graph.mat nodes edges seg offset

close( hwait )


