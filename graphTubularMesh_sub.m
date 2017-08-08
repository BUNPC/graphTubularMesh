

flagPlot = 0;
nIter = 0;
lstVmoveableHold = [];
while ~isempty( lstVmoveable )
    nIter = nIter + 1;
    
    % Display connected mesh
    if mod(nIter,1)==0 & flagPlot==1
        displayMeshTube( f,v, fConLeft, lstVcgHold, lstVmoveable, gV, seg, camera, 1, 3, lstVmoveableHold )
        pause(0.1)
    end
    
    
    camera.CP = campos;
    camera.CT = camtarget;
    camera.CVA = camva;
    camera.UP = camup;
    
    
    % get center point, remove end vertices and associated faces, and get new
    % end vertices and clean up
    
    % get center point
    nvg = nvg + 1;
    gV(nvg,:) = mean( v(lstVmoveable,:), 1 );
    
    % add to current graph segment
    seg(nSeg).nv = seg(nSeg).nv + 1;
    seg(nSeg).v( seg(nSeg).nv ) = nvg;
    dpos = [0 0 0];
    if seg(nSeg).nv > 1
        seg(nSeg).ne = seg(nSeg).ne + 1;
        seg(nSeg).e( seg(nSeg).ne, 1:2 ) = [ seg(nSeg).v( seg(nSeg).nv-1 ) seg(nSeg).v( seg(nSeg).nv ) ];
        dpos = gV(seg(nSeg).v( seg(nSeg).nv ),:) - gV(seg(nSeg).v( seg(nSeg).nv-1 ),:);
        dpos = dpos / norm(dpos);
    end

    % find end vertices that are propagating ahead of bisecting plane and
    % hold them from pruning this step
    if 0
        zv = zeros(length(lstVmoveable),1);
        for ii = 1:length(lstVmoveable)
            zv(ii) = norm( (v(lstVmoveable(ii),:) - gV(seg(nSeg).v( seg(nSeg).nv ),:) ) .* dpos );
        end
        lstVmoveableHold = lstVmoveable( find(zv>100) ); % instead of 2, use something like 2x median face edge length
        lstVmoveableCut = lstVmoveable( find(zv<=100) ); 
    end
    
    % remove end vertices
    vConLeft = setdiff( vConLeft, lstVmoveable );
%    vConLeft = setdiff( vConLeft, lstVmoveableCut );
    
    % remove corresponding faces
    for ii=1:length(lstVmoveable)
        [lst, ic] = find( f(fConLeft,:)==lstVmoveable(ii) );
%    for ii=1:length(lstVmoveableCut)
%        [lst, ic] = find( f(fConLeft,:)==lstVmoveableCut(ii) );
        fConLeft(lst) = [];
    end
    
    % new list of moveable vertices
    [ir,ic]=find(Medges(lstVmoveable,:)>0);
    lstVmoveable = intersect( setdiff(unique(ic),lstVmoveable), vConLeft);
%    [ir,ic]=find(Medges(lstVmoveableCut,:)>0);
%    ic = [ic; lstVmoveableHold];
%    lstVmoveable = intersect( setdiff(unique(ic),lstVmoveableCut), vConLeft);
    
    % prune any end vertices that don't have faces
    lstVmoveable = pruneVerticesWithNoFaces( lstVmoveable, f, fConLeft );
    
    % if reach end, then pick up from another start point and start a new
    % segment
    while isempty( lstVmoveable ) & ~isempty( lstVcgHold )
        lstVmoveable = lstVcgHold(end).v;
        
        lstVmoveable = pruneVerticesWithNoFaces( lstVmoveable, f, fConLeft );
        
        if seg(nSeg).nv>1
            nSeg = nSeg + 1;
        end
        seg(nSeg).v = lstVcgHold(end).nvg_con;
        seg(nSeg).nv = 1;
        seg(nSeg).e = [];
        seg(nSeg).ne = 0;

        lstVcgHold(end) = [];
    end
    
    % find connected groups
    % if more than 1, then start a new segment
    [nCG, lstVcg] = findConnectedGroups( lstVmoveable, Medges );
    if nCG==1
        lstVmoveable = lstVcg(1).v;
    elseif nCG>1
        nCG
        lstVmoveable = lstVcg(1).v;
        for ii=2:nCG
            lstVcgHold(end+1).v = lstVcg(ii).v;
            lstVcgHold(end).nvg_con = nvg;
        end
        
        if seg(nSeg).nv>1
            nSeg = nSeg + 1;
        end
        seg(nSeg).v = nvg;
        seg(nSeg).nv = 1;
        seg(nSeg).e = [];
        seg(nSeg).ne = 0;
    end

end

if seg(nSeg).nv==1
    seg(nSeg) = [];
    nSeg = nSeg - 1;
end
