function tex = snapGridToVertices(P, F, tex0, res, vertIdx)
    
    tex0 = res*matrixnormalize(tex0);
    
    nC = numel(vertIdx);
    if nC == 0
        tex = tex0;
        return
    end
    
    if iscolumn(vertIdx)
        vertIdx = vertIdx';
    end
    
    dim = size(P, 2);
    nP = size(P, 1);
    
    G = grad(P, F);
    
    % snap the first for now
    t = ceil(tex0(vertIdx(1), :)) - tex0(vertIdx(1), :);
    
    tex = tex0 + t;
    tex(vertIdx(2:end), :) = round(tex(vertIdx(2:end), :));
    
    
    roweq = 1:dim*nC;
    if (dim == 3)
        coleq = [vertIdx, nP+vertIdx, 2*nP+vertIdx];
        gradTex0 = [G*tex0(:, 1); G*tex0(:, 2); G*tex0(:, 3)];
        C = blkdiag(G, G, G);
    else
        coleq = [vertIdx, nP+vertIdx];
        gradTex0 = [G*tex0(:, 1); G*tex0(:, 2)];
        C = blkdiag(G, G);
    end
    
    
    idxeq = sub2ind([dim*nC, dim*nP], roweq, coleq);
    
    Aeq = zeros(dim*nC, dim*nP);
    Aeq(idxeq) = 1;
    
    beq = reshape(tex(vertIdx, :), [nC*dim 1]);
        
    tex = lsqlin(C, gradTex0,...
        [], [],...
        Aeq, beq,...
        zeros(dim*nP, 1), (1+res)*ones(dim*nP, 1), ...
        tex);
    
    tex = reshape(tex, [nP dim]);
end