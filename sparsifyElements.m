function dataOut = sparsifyElements(dataIn, stride, offset)
%% SPARSIFYELEMENTS Filters out elements to create a sparser truss graph.
%   Filters elements of input truss graph to only keep elements which trace
%   out isosurfaces with the value k*stride + offset, for non-negative
%   integers k. Both stride and offset can also be 3-tuples to specify a
%   different stride and offset value for each parameter. offset is
%   optional (defaults to zero).

    elem = dataIn.Elem;
    U = dataIn.NodeU;
    V = dataIn.NodeV;
    B = dataIn.NodeB;
    
    if nargin < 3
        offset = 0;
    end
    
    if numel(stride) == 1
        stride = [stride stride stride];
    end
    
    if numel(offset) == 1
        offset = [offset offset offset];
    end
    
    MOD = @(x, y) (x>0 & mod(x+offset(y), stride(y)) == 0);
    
    cU1 = MOD(U(elem(:, 1)), 1);
    cU2 = MOD(U(elem(:, 2)), 1);
    cV1 = MOD(V(elem(:, 1)), 2);
    cV2 = MOD(V(elem(:, 2)), 2);
    cB1 = B(elem(:, 1));
    cB2 = B(elem(:, 2));
    
    cU = (cU1 & U(elem(:, 1)) == U(elem(:, 2)));
    cV = (cV1 & V(elem(:, 1)) == V(elem(:, 2)));
    
    if size(dataIn.Node, 2) == 3
        W = dataIn.NodeW;
        cW1 = MOD(W(elem(:, 1)), 3);
        cW2 = MOD(W(elem(:, 2)), 3);
        cW = (cW1 & W(elem(:, 1)) == W(elem(:, 2)));
        
        cB = cB1 & cB2;
        cE = dataIn.ElemC;
        
        cUB = ((cU1 & cB2) | (cU2 & cB1)) & U(elem(:, 1)) == U(elem(:, 2));
        cVB = ((cV1 & cB2) | (cV2 & cB1)) & V(elem(:, 1)) == V(elem(:, 2));
        cWB = ((cW1 & cB2) | (cW2 & cB1)) & W(elem(:, 1)) == W(elem(:, 2));
        
        choiceBoundary = cE | ...   % feature edges
            (cU & (cVB | cB)) |...  % UV-lines interacting with the boundary
            (cV & (cWB | cB)) |...
            (cW & (cUB | cB));

        choiceInternal = (cU & cV) | (cV & cW) | (cW & cU);

        elem = elem(choiceBoundary | choiceInternal, :);
        
        if isfield(dataIn, 'Curves')
            curves = dataIn.Curves;
            idx = false(numel(curves), 1);
            for i=1:numel(curves)
                cE = curves(i).f;
                cU = MOD(curves(i).u, 1);
                cV = MOD(curves(i).v, 2);
                cW = MOD(curves(i).w, 3);
                cB = curves(i).b;
                
                idx(i) = (cB & (cU | cV | cW)) |... % Valid boundary curve: follows one of the desired isolines
                    ((cU & cV) | (cV & cW) | (cW & cU)) |... % Internal curve following 2 desired isolines
                    cE; % Feature curve
            end
            curves = curves(idx);
        end
    else
        cB = cB1 & cB2;
        elem = elem(cU | cV | cB, :);
    end
    
    dataOut = dataIn;
    dataOut.Elem = elem;
    if size(dataIn.Node, 2) == 3
        dataOut.ElemC = dataIn.ElemC(choiceBoundary | choiceInternal);
        if isfield(dataIn, 'Curves')
            dataOut.Curves = curves;
        end
    end
end