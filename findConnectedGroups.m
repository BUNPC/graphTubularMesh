function [nCG, lstVcg] = findConnectedGroups( lstVmoveable, Medges )

lstVmoveableLeft = lstVmoveable;
nCG = 0;
lstVcg = [];
lstVtmp = [];

while ~isempty(lstVmoveableLeft)
    
    nCG = nCG + 1;
    lstVcg(nCG).v = lstVmoveableLeft(1);
    
    ic = 1;
    while ~isempty(ic)
%        [ir,ic] = find(Medges(lstVcg(nCG).v(end),:)>0);
        ic = find(Medges(lstVcg(nCG).v(end),:)>0);
        ic = intersect(ic,lstVmoveableLeft);
        ic = setdiff(ic,lstVcg(nCG).v);
        if ~isempty(ic)
            lstVcg(nCG).v(end+1,1) = ic(1);
            if length(ic)>1
                lstVtmp(end+[1:(length(ic)-1)],1) = ic(2:end);
            end
        elseif ~isempty(lstVtmp)
            lstVtmp = unique(lstVtmp);
            ic = lstVtmp(end);
            if length(lstVtmp)>1
                lstVtmp(end) = [];
            else
                lstVtmp = [];
            end
            lstVcg(nCG).v(end+1,1) = ic(1);
        end
    end
    
    lstVmoveableLeft = setdiff( lstVmoveableLeft, lstVcg(nCG).v );

end