function lstVmoveable = pruneVerticesWithNoFaces( lstVmoveable, f, fConLeft )

    lstPrune = [];
    for ii=1:length(lstVmoveable)
        lst = find( f(fConLeft,:)==lstVmoveable(ii) );
        if isempty(lst)
            lstPrune(end+1) = ii;
        end
    end
    lstVmoveable(lstPrune) = [];
