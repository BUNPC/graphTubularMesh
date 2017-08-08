% To Do
% use imfill to mask out unconnected vessels
%   but what if node not connected to any vessels

function nodePos = centerNodesXYZ( nodePos, nodeEdges, I, centerStep1vox, nLst )

flagVisualize = 0;

if ~exist('centerStep1vox')
    centerStep1vox = 1;
end
if ~exist('nLst')
    nLst = [];
end
if isempty(nLst)
    nLst = 1:size(nodePos,1);
end

IThresh = 0.5;

nNodes = size(nodePos,1);
[ny,nx,nz] = size(I);


if ~flagVisualize
    hwait = waitbar(0,'Centering XYZ...');
end

nB = numNodeEdges( nodePos, nodeEdges );

nodeFlag = ones(size(nodePos,1),1);

nNodes_updated=length(nLst);
jNode=1;
nodePosTmp = nodePos;

while jNode<=nNodes_updated
    
    %jNode=1:length(nLst) %1:nNodes
    iNode = nLst(jNode);

    if ~flagVisualize
        waitbar(jNode/length(nLst),hwait);
    end

    if nB(iNode)<=2
        [lstE,lstC] = find(nodeEdges==iNode);

        if numel(lstE)==0,
            button = questdlg(['Failed at node ' num2str(iNode) ' (X,Y,Z) = ' num2str(nodePos(iNode,:)) ' Delete failing node?' ], 'Delete failing node?'); 
            if strcmp(button,'Yes'),
                nodeFlag(iNode) = 0;
%                [nodePos,im.nodeDiam,im.nodeDiamThetaIdx,im.nodeBC,im.nodeBCType,im.nodeType,im.nodeSegN,im.nodeEdges,im.edgeFlag] = removeNodes( nodeFlag, im.nodePos, im.nodeDiam, im.nodeDiamThetaIdx, im.nodeBC, im.nodeBCType, im.nodeType, im.nodeSegN,im.nodeEdges );                
                nNodes_updated=nNodes_updated-1;
%                continue; %LG not a good idea
                %[lstE,lstC] = find(im.nodeEdges==iNode); %LG I was getting
                %a bug with im.nodePos(iNode) at 47
            else
                keyboard;
            end;
        else
            pos1 = nodePos(iNode,:);
            pos2 = nodePos(nodeEdges(lstE(1),mod(lstC(1),2)+1),:);
            if length(lstE)==2
                pos0 = nodePos(nodeEdges(lstE(2),mod(lstC(2),2)+1),:);
            else
                pos0 = [];
            end
            
            r = norm(pos2-pos1);
            if r>0
                theta = acos((pos2(3)-pos1(3))/r);
                rho = norm(pos2(1:2)-pos1(1:2));
                phi = acos((pos2(1)-pos1(1))/rho);
                
                if~isempty(pos0)
                    r = norm(pos1-pos0);
                    if r>0
                        theta2 = acos((pos1(3)-pos0(3))/r);
                        rho = norm(pos1(1:2)-pos0(1:2));
                        phi2 = acos((pos1(1)-pos0(1))/rho);
                        
                        while (phi-phi2)>3.14159
                            phi2 = phi2 + 2*3.14159;
                        end
                        theta = (theta + theta2)/2;
                        phi = (phi+phi2)/2;
                    end
                end
                
                xLst = [-10:10];
                yLst = [-10:10];
                Isub = zeros(length(yLst),length(xLst));
                for ix = 1:length(xLst)
                    for iy = 1:length(yLst)
                        dpos = [xLst(ix) yLst(iy) 0]';
                        dpos = [[-cos(theta) 0 sin(theta)];[0 1 0];[-sin(theta) 0 cos(theta)]] * dpos;
                        dpos = [[cos(phi) -sin(phi) 0];[sin(phi) cos(phi) 0];[0 0 1]] * dpos;
                        iix = max(min(round(pos1(1)+dpos(1)),nx),1);
                        iiy = max(min(round(pos1(2)+dpos(2)),ny),1);
                        iiz = max(min(round(pos1(3)+dpos(3)),nz),1);
                        Isub(iy,ix) = I( iiy, iix, iiz );
                    end
                end
                
                Isub2 = Isub .* (Isub >= IThresh); 
                ic=11; ir=11; % Check if 11,11 is in vessel and correct if not
                r=1; theta=-30;
                while Isub2(ir,ic)==0 && r<5
                    theta = theta + 30; 
                    if theta>=360
                        theta = 0;
                        r=r+1;
                    end
                    ic = round(11+r*cos(theta*3.14159/180));
                    ir = round(11+r*sin(theta*3.14159/180));
                end
                if r>=5
                    ic=11; ir=11;
                end
                Isub = double(bwselect(Isub2==1,ic,ir,4)); % continuity condition here using bwselect. get rid of unconnected objects
                
                IsubSum = max(sum(Isub(:)),1)+eps;
                [xx,yy] = meshgrid(xLst,yLst);
                posM(1) = sum(xx(:).*Isub(:))/IsubSum;
                posM(2) = sum(yy(:).*Isub(:))/IsubSum;
                
                if flagVisualize
                    figure(1)
                    subplot(2,1,1)
                    imagesc(I(:,:,round(pos1(3))))
                    xlim([-20 20]+pos1(1));
                    ylim([-20 20]+pos1(2));
                    ht = text(pos1(1),pos1(2),'o');
                    subplot(2,1,2)
                    imagesc(xLst,yLst,Isub+Isub2,[0 2])
                    ht = text(ic-11,ir-11,'o');
                    ht = text(posM(1),posM(2),'x');
                    pause
                end
                
                
                if centerStep1vox
                    posM = posM / (max(norm(posM),1)+eps);
                end
                
                dpos = [posM(1) posM(2) 0]';
                dpos = [[-cos(theta) 0 sin(theta)];[0 1 0];[-sin(theta) 0 cos(theta)]] * dpos;
                dpos = [[cos(phi) -sin(phi) 0];[sin(phi) cos(phi) 0];[0 0 1]] * dpos;
                iix = max(min(pos1(1)+dpos(1),nx),1);
                iiy = max(min(pos1(2)+dpos(2),ny),1);
                iiz = max(min(pos1(3)+dpos(3),nz),1);
                nodePosTmp(iNode,:) = [iix iiy iiz];
                %        pause(0.1)
            end % end r>0
        end % end if numel(lstE)==0
    end  % end if nB(iNode)<=2
jNode=jNode+1;%LG to keep track when we remove nodes
end % end loop on iNode

lst = find(nodeFlag==1);
nodePos = nodePosTmp(lst,:);

close(hwait)
