function dataOut = collapseSmallEdges(dataIn, epsilon)
    v = dataIn.Node;
    e = dataIn.Elem;
    
    B = dataIn.NodeB;
    
    U = dataIn.NodeU;
    V = dataIn.NodeV;
    
    dim = size(v, 2);
    
    if dim == 3
        EC = dataIn.Elem(dataIn.ElemC, :);
        F = accumarray(EC(:), 1) > 0;
        F = [F; false(numel(B) - numel(F), 1)];
        W = dataIn.NodeW;
        EC = dataIn.ElemC;
        deg = accumarray([e(EC, 1); e(EC, 2)], 1);
        deg = [deg; zeros(size(v, 1) - numel(deg), 1)];
%         C = dataIn.NodeC;
        C = deg > 0;
    else
        C = false(size(B));
        F = C;
    end
    
    eLen = sqrt( sum( (v(e(:, 1), :) - v(e(:, 2), :)).^2, 2 ) );
    
    % very slow implementation: removes the (one) smallest edge in every 
    % iteration
    while min(eLen) < epsilon
        [~, idx] = min(eLen);
        pt = e(idx, :);
        e(idx, :) = [];
        if dim==3
            EC(idx, :) = [];
        end
        
        % "Priority" of vertices to keep their input positions: First,
        % feature points. Second, boundary points. Then, the rest.
        
        % if one is a corner, and the other is not
        if C(pt(1)) && ~C(pt(2))
            v(pt(2), :) = v(pt(1), :);
            U(pt(1)) = max(U(pt));
            V(pt(1)) = max(V(pt));
            B(pt(1)) = max(B(pt));
            if dim == 3
                C(pt(1)) = max(C(pt));
                W(pt(1)) = max(W(pt));
            end
            e(e==pt(2)) = pt(1);
        elseif C(pt(2)) && ~C(pt(1))
            v(pt(1), :) = v(pt(2), :);
            U(pt(2)) = max(U(pt));
            V(pt(2)) = max(V(pt));
            B(pt(2)) = max(B(pt));
            if dim == 3
                C(pt(2)) = max(C(pt));
                W(pt(2)) = max(W(pt));
            end
            e(e==pt(1)) = pt(2);
        % if one lies on a feature curve, and the other does not
        elseif F(pt(1)) && ~F(pt(2))
            v(pt(2), :) = v(pt(1), :);
            U(pt(1)) = max(U(pt));
            V(pt(1)) = max(V(pt));
            B(pt(1)) = max(B(pt));
            if dim == 3
                C(pt(1)) = max(C(pt));
                W(pt(1)) = max(W(pt));
            end
            e(e==pt(2)) = pt(1);
        elseif F(pt(2)) && ~F(pt(1))
            v(pt(1), :) = v(pt(2), :);
            U(pt(2)) = max(U(pt));
            V(pt(2)) = max(V(pt));
            B(pt(2)) = max(B(pt));
            if dim == 3
                C(pt(2)) = max(C(pt));
                W(pt(2)) = max(W(pt));
            end
            e(e==pt(1)) = pt(2);
        % if one is a boundary point and the other is not
        elseif B(pt(1)) && ~B(pt(2))
            v(pt(2), :) = v(pt(1), :);
            U(pt(1)) = max(U(pt));
            V(pt(1)) = max(V(pt));
            B(pt(1)) = max(B(pt));
            if dim == 3
                C(pt(1)) = max(C(pt));
                W(pt(1)) = max(W(pt));
            end
            e(e==pt(2)) = pt(1);
        elseif B(pt(2)) && ~B(pt(1))
            v(pt(1), :) = v(pt(2), :);
            U(pt(2)) = max(U(pt));
            V(pt(2)) = max(V(pt));
            B(pt(2)) = max(B(pt));
            if dim == 3
                C(pt(2)) = max(C(pt));
                W(pt(2)) = max(W(pt));
            end
            e(e==pt(1)) = pt(2);
        % if both points have the same "priority", take the average
        % position
        else
            v(pt(1), :) = (v(pt(1), :) + v(pt(2), :))/2;
            U(pt(2)) = max(U(pt));
            V(pt(2)) = max(V(pt));
            B(pt(2)) = max(B(pt));
            if dim == 3
                C(pt(2)) = max(C(pt));
                W(pt(2)) = max(W(pt));
            end
            e(e==pt(1)) = pt(2);
        end
             
        eLen = sqrt( sum( (v(e(:, 1), :) - v(e(:, 2), :)).^2, 2 ) );
    end
    
    [v, I, J] = remove_unreferenced(v, e);
    e = I(e);
    
    B = B(J);
    U = U(J);
    V = V(J);
    
    if dim==3
        C = C(J);
        W = W(J);
    end

    e = sort(e, 2);
    [e, I] = unique(e, 'rows');
    
    dataOut = dataIn;
    dataOut.Node = v;
    dataOut.Elem = e;
    dataOut.NodeB = B;
    dataOut.NodeU = U;
    dataOut.NodeV = V;
    
    if dim==3
        EC = EC(I);
        dataOut.ElemC = EC;
        dataOut.NodeC = C;
        dataOut.NodeW = W;
    end
end