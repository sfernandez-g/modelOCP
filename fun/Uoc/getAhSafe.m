function Ah = getAhSafe(cap, ind)
    if isempty(ind)
        Ah = 0;
    else
        Ah = cap(ind);
        Ah = Ah - Ah(1);
    end