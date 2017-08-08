function nB = numNodeEdges( nodePos, nodeEdges )

nNodes = size(nodePos,1);
nB=zeros(nNodes,1);
hwait = waitbar( 0,'Calculating Number Edges Per Node' );
for ii=1:nNodes
    nB(ii)=length(find(nodeEdges(:,1)==ii | nodeEdges(:,2)==ii));
end
close(hwait);


