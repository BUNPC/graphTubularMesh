function nodePos = centerNodes( nodePos, nodeEdges, I )
% I is y,x,z and nodePos is x,y,z

% im.nX; im.nY; im.nZ;
% im.Hvox(1:3)
% im.I - use segmented image on which mesh is based
% Ithresh

flagVisualize = 0;
Hvox = [1 1 1]; % voxel size ?
Ithresh = 0.5;
nLst = 1:size(nodePos,1);
centerStep = 1;

% ch=1, 2, 3 from 'XY','XY & estimate diameter','Z'
ch = 1;

% Center Nodes in XY or Z
minDiam = 5;
maxDiam = 90;


nN = size(nodePos,1);

if ch==1 | ch==2
    pIdx = [1 2]; % INDICES OF X AND Y in I
    nnX = size(I,2); %im.nX;
    nnY = size(I,1); %im.nY;
elseif ch==3
    pIdx = [2 3]; % INDICES OF Z AND Y in I
    nnX = size(I,3); %im.nZ;
    nnY = size(I,1); %im.nY;
end

nx = ceil(maxDiam / Hvox(pIdx(1)));
ny = ceil(maxDiam / Hvox(pIdx(2)));
%nz = ceil(maxDiam / im.Hvox(3));
hx = Hvox(pIdx(1));
hy = Hvox(pIdx(2));
%hz = im.Hvox(3);


[nr,nc,ns] = size(I);
if ch==1 | ch==2
    [nr,nc] = size(I(:,:,1));
elseif ch==3
    [nr,nc] = size(squeeze(I(:,1,:)));
end

% create mapping of different angular cuts
iTheta = 0;
if ch==1 | ch==2
    thetaLst = [0:pi/10:pi-.1];
elseif ch==3
    thetaLst = 0;
end
nTheta = length(thetaLst);
rhoLst = [-maxDiam:1:maxDiam];
nRho = length(rhoLst);
lineMap = zeros(nRho,nTheta);
for iTheta = 1:nTheta
    for iRho = 1:nRho
        xx = rhoLst(iRho) * cos(thetaLst(iTheta));
        yy = rhoLst(iRho) * sin(thetaLst(iTheta));
        
        ix = round((xx+maxDiam)/hx);
        iy = round((yy+maxDiam)/hy);
        
        lineMap(iRho,iTheta) = (2*ny+1)*ix + iy + 1;
    end
end


% Loop over nodes
if ~flagVisualize
    hwait = waitbar(0,'Centering...');
end
for iii=1:length(nLst) %1:size(im.nodePos,1)
    ii = nLst(iii);
    if ~flagVisualize
        waitbar(iii/length(nLst),hwait);
    end
    
    eLst = find(nodeEdges(:,1)==ii | nodeEdges(:,2)==ii);
    
    % if centering on Z, verify that mean edge direction is not in z
    if ch==3
        eLst = find(nodeEdges(:,1)==ii | nodeEdges(:,2)==ii);
        c3 = 0;
        for jj=1:length(eLst)
            pos1 = nodePos(nodeEdges(eLst(jj),1),:);
            pos2 = nodePos(nodeEdges(eLst(jj),2),:);
            c3 = c3 + abs((pos1(3)-pos2(3)) / norm(pos1-pos2));
        end
        c3 = c3 / length(eLst);
    else
        c3 = 0;
    end
    
    % if XY centering or (z centering and c3<Thresh) then center
    if c3<0.5
        pos=round(nodePos(ii,:));
        
        if (pos(pIdx(1))-nx)>=1
            dx1 = 1;
            sx1 = pos(pIdx(1))-nx;
        else
            dx1 = nx-pos(pIdx(1))+2;
            sx1 = 1;
        end
        if (pos(pIdx(1))+nx)<=nc
            dx2 = 2*nx+1;
            sx2 = pos(pIdx(1))+nx;
        else
            dx2 = 2*nx+1-((pos(pIdx(1))+nx)-nc);
            sx2 = nc;
        end
        
        if (pos(pIdx(2))-ny)>=1
            dy1 = 1;
            sy1 = pos(pIdx(2))-ny;
        else
            dy1 = ny-pos(pIdx(2))+2;
            sy1 = 1;
        end
        if (pos(pIdx(2))+ny)<=nr
            dy2 = 2*ny+1;
            sy2 = pos(pIdx(2))+ny;
        else
            dy2 = 2*ny+1-((pos(pIdx(2))+ny)-nr);
            sy2 = nr;
        end
        
        Id = 32*ones(2*ny+1, 2*nx+1);
        if ch==1 | ch==2
            Id(dy1:dy2,dx1:dx2) = I(sy1:sy2,sx1:sx2,min(pos(3),ns));
        elseif ch==3
            Id(dy1:dy2,dx1:dx2) = squeeze(I(sy1:sy2,pos(1),sx1:sx2));
        end
        Id(find(Id<Ithresh))=0;
        
        
        I1 = Id(lineMap)<Ithresh;
        I2 = [];
        for ii2=1:nTheta
            I2(:,ii2) = imfill(I1(:,ii2), round(nRho/2), [0 1 0;0 1 0;0 1 0]) - I1(:,ii2);
        end
        
        [diam,iTheta] = min(sum(I2,1));
        % If diam = 0 the try to extrapolate out a few voxels to find the
        % vessel
        if diam == 0
            for iOff = 1:3  % this will only extrapolate +/- 3 voxels
                for ii2=1:nTheta
                    I2(:,ii2) = imfill(I1(:,ii2), round(nRho/2)+iOff, [0 1 0;0 1 0;0 1 0]) - I1(:,ii2);
                end
                lstGT0 = find(sum(I2,1)>0);
                if ~isempty(lstGT0)
                    [diam,iThetaTmp] = min(sum(I2(:,lstGT0),1));
                    iTheta = lstGT0(iThetaTmp);
                    break
                end
                for ii2=1:nTheta
                    I2(:,ii2) = imfill(I1(:,ii2), round(nRho/2)-iOff, [0 1 0;0 1 0;0 1 0]) - I1(:,ii2);
                end
                lstGT0 = find(sum(I2,1)>0);
                if ~isempty(lstGT0)
                    [diam,iThetaTmp] = min(sum(I2(:,lstGT0),1));
                    iTheta = lstGT0(iThetaTmp);
                    break
                end
            end
        end
        if ch==3
            iTheta = 1;
        end
        theta = thetaLst(iTheta);
        lstDiam = find(I2(:,iTheta)==1);
        
        if ~isempty(lstDiam)
            dR = mean(lstDiam) - maxDiam;
        else
            dR = 0;
        end
        if centerStep==1
            dR = sign(dR)*min(abs(dR),1);
        end
        %%%%%%%%%%%%%%%
        % check if this moves node off of 1:nx or 1:ny
        % if so then set dR = 0
        foo = (nodePos(ii,pIdx(1))+cos(theta)*dR/hx);
        if foo<1 | foo>nnX
            dR = 0;
        end
        foo = (nodePos(ii,pIdx(2))+sin(theta)*dR/hy);
        if foo<1 | foo>nnY
            dR = 0;
        end
        
        if mod(ii,10)==1 & flagVisualize
            if gcf~=12
                figure(12);
            end
            subplot(1,3,1);
            if ch==1 | ch==2
                imagesc(I(:,:,min(pos(3),ns)));
            elseif ch==3
                imagesc(squeeze(I(:,pos(1),:)));
            end
            axis image
            ht=text(double(pos(pIdx(1))),double(pos(pIdx(2))),'x');
            set(ht,'color','m');
            set(ht,'fontweight','bold');
            set(ht,'horizontalalignment','center')
            
            Id2 = Id;
            Id2(lineMap(lstDiam,iTheta)) = 0;
            subplot(1,3,2)
            imagesc(Id2,[0 32])
            axis image
            ht=text(nx+1,ny+1,'x');
            set(ht,'color','k');
            set(ht,'fontweight','bold');
            set(ht,'horizontalalignment','center')
            ht=text(nx+1+cos(theta)*dR/hx,ny+1+sin(theta)*dR/hy,'x');
            set(ht,'color','m');
            set(ht,'fontweight','bold');
            set(ht,'horizontalalignment','center')
            title( sprintf('Node %d',ii) )
            
            cm = jet(32);
            cm(1,:) = [1 1 1];
            colormap(cm)
            
            subplot(1,3,3)
            imagesc(I2+(1-I1))
            
            drawnow
            pause
        end
        
        if ch==2
            im.nodeDiamEst(ii) = diam;
            im.nodeDiam(ii) = diam;
            im.nodeDiamThetaIdx(ii) = iTheta;
        end
        nodePos(ii,pIdx) = nodePos(ii,pIdx) + [cos(theta)/hx sin(theta)/hy]*dR;
    end % end check on c3 thresh for z centering
end % end loop over nodes
if ~flagVisualize
    close(hwait);
end

