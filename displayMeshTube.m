function displayMeshTube( f,v, fConLeft, lstVcgHold, lstVmoveable, gV, seg, camera, edgeAlpha, fig, lstVmoveableHold )

col = 'mcyrg';
colSeg = 'mcbrg';
nCGhold = length(lstVcgHold);

if ~exist('fig')
    fig = 3;
end
if isempty(fig)
    fig = 3;
end

figure(fig);

% display the mesh
h=trisurf(f(fConLeft,:),v(:,1),v(:,2),v(:,3),'facecolor','none','edgecolor','b');
set(h,'EdgeAlpha',edgeAlpha)
hold on

% display vertex groups on hold
if nCGhold>0
    for ii=1:nCGhold
        h = plot3( v(lstVcgHold(ii).v,1), v(lstVcgHold(ii).v,2), v(lstVcgHold(ii).v,3), [col(mod(ii-1,4)+2) 'x'] );
        set(h,'markersize',6)
    end
end

% display the moveable vertex group
h=plot3(v(lstVmoveable,1), v(lstVmoveable,2), v(lstVmoveable,3), 'm.');
set(h,'markersize',20)
if exist('lstVmoveableHold')
    h=plot3(v(lstVmoveableHold,1), v(lstVmoveableHold,2), v(lstVmoveableHold,3), 'g.');
    set(h,'markersize',20)
end

% display the graph
% if ~isempty(gV)
%     h=plot3(gV(:,1), gV(:,2), gV(:,3), 'b.');
%     set(h,'markersize',20)
% end

% display the graph segments
for ii = 1:length(seg)
    if ~isempty( seg(ii).v )
        h = plot3( gV(seg(ii).v,1), gV(seg(ii).v,2), gV(seg(ii).v,3), [colSeg(mod(ii-1,5)+1) '.-'] );
        set(h,'markersize',20)
    end
end

hold off

% set the camera and lighting
daspect([1,1,1])
%view(3); 
axis tight
camlight
lighting gouraud

if ~isempty(camera.CP)
    campos(camera.CP);
    camtarget(camera.CT);
    camva(camera.CVA);
    camup(camera.UP);
end