function dataOut = tex2CurvesTet(dataIn, res, boundaryOff, snapVertIdx)
% Trace the isolines of a texture parametrization to generate a truss
% layout
%
% Input
% dataIn: struct with the mesh description (dataIn.V, dataIn.T) and the 
% parametrization (dataIn.u, dataIn.v, dataIn.w)
% res: the resolution parameter (`rho` in the paper)
% boundaryOff: If true, then does not trace boundary elements
% snapVertIdx: Optional. Set of indices into dataIn.V to snap the truss
% layout vertices to. That is, a truss node will lie on each of the indexed
% vertices


% Single-character abbreviations:
% (U, V, W): the three texture coordinates
% F: Faces (triangles)
% T: Simplices (tris for 2D, tets for 3D)
% E: Edges
% P: Vertices/positions/points
% I: integer
% (X, Y, Z): the three world coordinates
% B: Domain boundary

    if nargin < 3
        boundaryOff = false;
    end
    if nargin < 4
        snapVertIdx = [];
    end
    
    DONOTTRACEINTERNAL = false; % removes both nodes and elements
    DONOTTRACEBOUNDARY = boundaryOff;  % removes elements only
	
	
    P = dataIn.V;
    T = sort(dataIn.T, 2);
    % get F from T by enumerating all 4 choices on each tet, and then
    % taking a set union by removing duplicates
    choice = nchoosek(1:4, 3);
    F = unique(sort([...
        T(:, choice(1, :));...
        T(:, choice(2, :));...
        T(:, choice(3, :));...
        T(:, choice(4, :))...
        ], 2), 'rows');
    
	% Get E from T
    E = edges(T);
    
    nP = size(P, 1);
    nT = size(T, 1);
    nF = size(F, 1);
    nE = size(E, 1);
    
    tex = dataIn.u;
    
    adjacency = adjacency_list(E);
    
    tex = 1 + snapGridToVertices(P, T, tex, res, snapVertIdx);
    
    epsLarge = 1e-7;
    epsSmall = 1e-10;
    
    oneRingMax = [cellfun(@(x) max(tex(x, 1)), adjacency),...
        cellfun(@(x) max(tex(x, 2)), adjacency), ...
        cellfun(@(x) max(tex(x, 3)), adjacency)];
    texRound = round(tex);
    idxTexInt = abs(texRound-tex) < epsSmall;
    idxTexLocalMax = tex > oneRingMax - epsLarge;
    tex(idxTexInt & ~idxTexLocalMax) = ...
        tex(idxTexInt & ~idxTexLocalMax) - epsLarge;
    tex(idxTexInt & idxTexLocalMax) = ...
        tex(idxTexInt & idxTexLocalMax) + epsLarge;

    U = tex(:, 1);
    V = tex(:, 2);
    W = tex(:, 3);
    
    clear tex;
    
    
    [~, E_F1] = ismember(F(:, [1, 2]), E, 'rows');
    [~, E_F2] = ismember(F(:, [2, 3]), E, 'rows');
    [~, E_F3] = ismember(F(:, [1, 3]), E, 'rows');
    
    E_F = [E_F1, E_F2, E_F3];
        
    [~, F_T1] = ismember(T(:, [1 2 3]), F, 'rows');
    [~, F_T2] = ismember(T(:, [1 2 4]), F, 'rows');
    [~, F_T4] = ismember(T(:, [1 3 4]), F, 'rows');
    [~, F_T3] = ismember(T(:, [2 3 4]), F, 'rows');

    F_T = [F_T1, F_T2, F_T3, F_T4];

    isFBoundary = is_boundary_facet(F, T);
    
    F_E = cellfun(...
        @(x) find(isFBoundary & (E_F1 == x | E_F2 == x | E_F3 == x)),...
        num2cell((1:nE)'),...
        'UniformOutput', false);
    
    clear E_F1 E_F2 E_F3;
    clear F_T1 F_T2 F_T3 F_T4;
    
    % Points of intersection between an edge and an integral isosurface
    PUI_E = cell(nE, 1);
    PVI_E = PUI_E;
    PWI_E = PUI_E;
    
    % Values of the possibly non-integral tex-coordinates on these points
    VonPUI_E = PUI_E;
    WonPUI_E = PUI_E;
    UonPVI_E = PUI_E;
    WonPVI_E = PUI_E;
    UonPWI_E = PUI_E;
    VonPWI_E = PUI_E;
    
    % Values of the integral tex-coordinates on these points
    UI_E = cell(nE, 1);
    VI_E = UI_E;
    WI_E = UI_E;
    
    % Points of intersection between integral isocurves (curves with two
    % integral tex-coordiantes) and faces (tris)
    PUVI_F = cell(nF, 1);
    PVWI_F = PUVI_F;
    PWUI_F = PUVI_F;
    
    % Values of the third tex-coordinate on these points
    WonPUVI_F = PUVI_F;
    UonPVWI_F = PUVI_F;
    VonPWUI_F = PUVI_F;
    
    % Values of the integral tex-coordinates on these points
    U_UVI_F = PUVI_F;
    V_UVI_F = PUVI_F;
    V_VWI_F = PUVI_F;
    W_VWI_F = PUVI_F;
    W_WUI_F = PUVI_F;
    U_WUI_F = PUVI_F;
    
    allocSize = 2*1e6;
    % For final output: nodes (graph vertices)
    node = zeros(allocSize, 3);
    nodeSize = 0;
    % node properties
    nodeIsBoundary = false(allocSize, 1);
    nodeUI = int32(-ones(allocSize, 1));
    nodeVI = nodeUI;
    nodeWI = nodeUI;
    
    nodeU = -inf(allocSize, 1);
    nodeV = nodeU;
    nodeW = nodeU;
    
    nodeSizeTillPrevFace = zeros(nF, 1);
    nodeSizeTillPrevTet = zeros(nT, 1);
    % stores the edge idx for nodes lying on edges
    E_node = int32(zeros(allocSize, 1));
    % face idx for nodes lying on faces
    F_node = E_node;
    T_node = E_node;
    
    isEBoundary = false(nE, 1);
    for i=1:3
        isEBoundary(E_F(isFBoundary, i)) = true;
    end
    
    % For final output: elements (graph edges)
    elem = zeros(2*allocSize, 2);
    elemSize = 0;
    
    boundaryAngleCosine_E = ones(nE, 1);
    for i = 1:nE
        if isEBoundary(i)
            a = E(i, 1);
            b = E(i, 2);
            edge = P(a, :) - P(b, :);
            faces = F_E{i};
            
            c1 = setdiff(F(faces(1), :), E(i, :));
            c2 = setdiff(F(faces(2), :), E(i, :));
            n1 = cross(edge, P(c1, :) - P(b, :));
            n2 = cross(P(c2, :) - P(b, :), edge);
            boundaryAngleCosine_E(i) = dot(n1, n2)/(norm(n1) * norm(n2));
        end
    end
    
    isEBCurved = boundaryAngleCosine_E < 0.9;   % is E on a curved part of the boundary

    disp(['# of boundary edges: ' num2str(sum(isEBCurved))]);
    
    
    isPBoundary = false(nP, 1);
    isPBoundary(E(isEBoundary, 1)) = true;
    isPBoundary(E(isEBoundary, 2)) = true;
    
    EBCurved = E(isEBCurved, :);
    numIncidentEBCurved_P = accumarray(EBCurved(:), 1);
    numIncidentEBCurved_P = [numIncidentEBCurved_P; zeros(nP - max(EBCurved(:)), 1)];
    
    boundaryAngleCosine_P = ones(nP, 1);
    for i = 1:nP
        if isPBoundary(i) && numIncidentEBCurved_P(i) >= 2
            neighbours = adjacency{i}(isPBoundary(adjacency{i}));
            n1 = -1;
            n2 = -1;
            for n = neighbours
                v = sort([i, n]);
                % find the sharp edge {i, n}
                if numel(find(E(:, 1)==v(1) & E(:, 2)==v(2) & isEBCurved))
                    if (n1 == -1)
                        n1 = n;
                    else
                        n2 = n;
                        break;
                    end
                end
            end
            vec1 = P(i, :) - P(n1, :);
            vec2 = P(n2, :) - P(i, :);
            boundaryAngleCosine_P(i) = dot(vec1, vec2)/(norm(vec1) * norm(vec2));
        end
    end
    
    isPBCurved = boundaryAngleCosine_P < 0.5;   % is P on a highly curved part of the boundary
    numPBCurved = sum(isPBCurved);
    disp(['# of boundary corners: ' num2str(numPBCurved)]);
    
    
    for i = 1:nE
%         disp(['Edge # ' num2str(i)]);
        edge = E(i, :);
        
        % find points on intersections b/w integral u-isosurfaces and edges
        [U12, p1Idx] = sort([U(edge(1)), U(edge(2))]);
        edge = edge(p1Idx);
        
        p1 = P(edge(1), :);
        p2 = P(edge(2), :);
        
        v1 = V(edge(1));
        v2 = V(edge(2));
        
        w1 = W(edge(1));
        w2 = W(edge(2));
        
        UI_E{i} = integersOpen(U12(1), U12(2));
        
        alphaU = (UI_E{i} - U12(1)) / (U12(2) - U12(1));
        

        PUI_E{i} = (1 - alphaU)*p1 + alphaU*p2;
        VonPUI_E{i} = (1 - alphaU)*v1 + alphaU*v2;
        WonPUI_E{i} = (1 - alphaU)*w1 + alphaU*w2;
        
        % if edge is on the boundary, add these points to nodes
        if isEBoundary(i)
            n = numel(alphaU);
            nodeIdx = nodeSize + (1:n);
            node(nodeIdx, :) = PUI_E{i};
            nodeUI(nodeIdx) = UI_E{i};
            nodeU(nodeIdx) = UI_E{i};
            nodeV(nodeIdx) = VonPUI_E{i};
            nodeW(nodeIdx) = WonPUI_E{i};
            E_node(nodeIdx) = i;
            nodeIsBoundary(nodeIdx) = true;
            nodeSize = nodeSize + n;
        end
        
        % do the same for integral v-isosurfaces
        [V12, p1Idx] = sort([V(edge(1)), V(edge(2))]);
        edge = edge(p1Idx);
        
        p1 = P(edge(1), :);
        p2 = P(edge(2), :);
        
        u1 = U(edge(1));
        u2 = U(edge(2));
        
        w1 = W(edge(1));
        w2 = W(edge(2));
        
        VI_E{i} = integersOpen(V12(1), V12(2));
        
        alphaV = (VI_E{i} - V12(1)) / (V12(2) - V12(1));
        
        
        PVI_E{i} = (1 - alphaV)*p1 + alphaV*p2;
        UonPVI_E{i} = (1 - alphaV)*u1 + alphaV*u2;
        WonPVI_E{i} = (1 - alphaV)*w1 + alphaV*w2;
        
        if isEBoundary(i)
            n = numel(alphaV);
            nodeIdx = nodeSize + (1:n);
            node(nodeIdx, :) = PVI_E{i};
            nodeVI(nodeIdx) = VI_E{i};
            nodeU(nodeIdx) = UonPVI_E{i};
            nodeV(nodeIdx) = VI_E{i};
            nodeW(nodeIdx) = WonPVI_E{i};
            E_node(nodeIdx) = i;
            nodeIsBoundary(nodeIdx) = true;
            nodeSize = nodeSize + n;
        end
        
        % finally, do this for integral w-isosurfaces
        [W12, p1Idx] = sort([W(edge(1)), W(edge(2))]);
        edge = edge(p1Idx);
        
        p1 = P(edge(1), :);
        p2 = P(edge(2), :);
        
        u1 = U(edge(1));
        u2 = U(edge(2));
        
        v1 = V(edge(1));
        v2 = V(edge(2));
        
        WI_E{i} = integersOpen(W12(1), W12(2));
        
        alphaW = (WI_E{i} - W12(1)) / (W12(2) - W12(1));
        
        
        PWI_E{i} = (1 - alphaW)*p1 + alphaW*p2;
        UonPWI_E{i} = (1 - alphaW)*u1 + alphaW*u2;
        VonPWI_E{i} = (1 - alphaW)*v1 + alphaW*v2;
        
        if isEBoundary(i)
            n = numel(alphaW);
            nodeIdx = nodeSize + (1:n);
            node(nodeIdx, :) = PWI_E{i};
            nodeWI(nodeIdx) = WI_E{i};
            nodeU(nodeIdx) = UonPWI_E{i};
            nodeV(nodeIdx) = VonPWI_E{i};
            nodeW(nodeIdx) = WI_E{i};
            E_node(nodeIdx) = i;
            nodeIsBoundary(nodeIdx) = true;
            nodeSize = nodeSize + n;
        end
    end
    
%     disp('# boundary edges:');
%     disp(elemSize);
    
    % State of final graph so far:
    % node-> all nodes lying on boundary edges
    % elem-> empty
    
    for i=1:nF
        nodeSizeTillPrevFace(i) = nodeSize;
        % look at all U-lines on all faces, and find intersection points with
        % V-lines and with W-lines. For boundary faces, also build the
        % connectivity graph by adding elements forming the U-lines
%         disp(['Face # ' num2str(i)]);
%         disp('U-lines');
        for j=1:3
            e1 = E_F(i, mod(j+1, 3)+1);
            e2 = E_F(i, mod(j+2, 3)+1);
            
            minUI = max(min(UI_E{e1}), min(UI_E{e2}));
            maxUI = min(max(UI_E{e1}), max(UI_E{e2}));
            
            % idx into UIonE{e1} and UIonE{e2}
            [~, shift1] = ismember(minUI, UI_E{e1});
            if isempty(shift1)
                shift1 = 0;
            end
            idx1 = shift1 + (0:maxUI-minUI);
            [~, shift2] = ismember(minUI, UI_E{e2});
            if isempty(shift2)
                shift2 = 0;
            end
            idx2 = shift2 + (0:maxUI-minUI);
            
            for k=1:numel(idx1)
                % We have to interpolate b/w the vertices at idx1(k)
                % position in PUI_E{e1} and idx2(k) in PUI_E{e2}.
                % UI_E{e1/2} and VonUI_E{e1/2} follow the same indexing
                k1 = idx1(k);
                k2 = idx2(k);
                
                ePair = [e1 e2];
                
                P12 = [PUI_E{e1}(k1, :); PUI_E{e2}(k2, :)];
                
                u = UI_E{e1}(k1);   % U is the same for both

                V12 = [VonPUI_E{e1}(k1); VonPUI_E{e2}(k2)];
                W12 = [WonPUI_E{e1}(k1); WonPUI_E{e2}(k2)];
                
                % point with smaller v or w is at alpha = 0
                if abs(diff(V12)) > abs(diff(W12))
                    [~, idxSort] = sort(V12);
                    sorter = 2;
                else
                    [~, idxSort] = sort(W12);
                    sorter = 3;
                end
                V12 = V12(idxSort);
                P12 = P12(idxSort, :);
                ePair = ePair(idxSort);
                W12 = W12(idxSort);
                
                if V12(1) > V12(2)
                    rangeV = integersHalfOpen(V12(2), V12(1));
                    alphaV = (rangeV - V12(2)) / (V12(1) - V12(2));
                    alphaV = 1 - alphaV;
                else
                    rangeV = integersHalfOpen(V12(1), V12(2));
                    alphaV = (rangeV - V12(1)) / diff(V12);
                end
                
                if W12(1) > W12(2)
                    rangeW = integersHalfOpen(W12(2), W12(1));
                    alphaW = (rangeW - W12(2)) / (W12(1) - W12(2));
                    alphaW = 1 - alphaW;
                else
                    rangeW = integersHalfOpen(W12(1), W12(2));
                    alphaW = (rangeW - W12(1)) / diff(W12);
                end
                
                nvi = numel(alphaV);
                viNodes = (1 - alphaV)*P12(1, :) + alphaV*P12(2, :);
                w_vi = (1 - alphaV)*W12(1) + alphaV*W12(2);
                PUVI_F{i} = [PUVI_F{i}; viNodes];
                U_UVI_F{i} = [U_UVI_F{i}; repmat(u, nvi, 1)];
                V_UVI_F{i} = [V_UVI_F{i}; rangeV];
                WonPUVI_F{i} = [WonPUVI_F{i}; w_vi];

                nodeIdx_v = nodeSize + (1:nvi)';
                node(nodeIdx_v, :) = viNodes;
                nodeUI(nodeIdx_v) = u;
                nodeVI(nodeIdx_v) = rangeV;
                nodeU(nodeIdx_v) = u;
                nodeV(nodeIdx_v) = rangeV;
                nodeW(nodeIdx_v) = w_vi;
                F_node(nodeIdx_v) = i;
                
                nwi = numel(alphaW);
                wiNodes = (1 - alphaW)*P12(1, :) + alphaW*P12(2, :);
                v_wi = (1 - alphaW)*V12(1) + alphaW*V12(2);
                PWUI_F{i} = [PWUI_F{i}; wiNodes];
                U_WUI_F{i} = [U_WUI_F{i}; repmat(u, nwi, 1)];
                W_WUI_F{i} = [W_WUI_F{i}; rangeW];
                VonPWUI_F{i} = [VonPWUI_F{i}; v_wi];

                nodeIdx_w = nodeSize + nvi + (1:nwi)';
                node(nodeIdx_w, :) = wiNodes;
                nodeUI(nodeIdx_w) = u;
                nodeWI(nodeIdx_w) = rangeW;
                nodeU(nodeIdx_w) = u;
                nodeV(nodeIdx_w) = v_wi;
                nodeW(nodeIdx_w) = rangeW;
                F_node(nodeIdx_w) = i;
                
                % add elements on the boundary
                if isFBoundary(i)
                    % find existing nodes on edges
                    e1NodeIdx = find(E_node==ePair(1) & nodeUI==u);
                    e2NodeIdx = find(E_node==ePair(2) & nodeUI==u);

                    % combine vars for vi and wi and sort by v
                    v = [double(rangeV); v_wi];
                    w = [w_vi; double(rangeW)];
                    nodeIdx = [nodeIdx_v; nodeIdx_w];
                    n = nvi + nwi;
                    if sorter==2
                        [~, sortIdx] = sort(v);
                    else
                        [~, sortIdx] = sort(w);
                    end
                    nodeIdx = nodeIdx(sortIdx);
                    
                    nodeIsBoundary(nodeIdx) = true;

                    
                    if n
                        % add the internal elems
                        elem(elemSize + (1:n-1), :) = [nodeIdx(1:n-1), nodeIdx(2:n)];
                        elemSize = elemSize + n-1;
                        
                        % add two elems if there is any internal node
                        elem(elemSize + (1:2), :) = [...
                            e1NodeIdx, nodeIdx(1);
                            nodeIdx(end), e2NodeIdx];
                        elemSize = elemSize + 2;
                    else
                        % otherwise, just add the elem connecting the two boundary
                        % nodes
                        elem(elemSize + 1, :) = [e1NodeIdx, e2NodeIdx];
                        elemSize = elemSize + 1;
                    end
                end
                
                nodeSize = nodeSize + nvi + nwi;
            end    
        end
        
        % Repeat for V-lines now. But we now have to search for UV points
        % in the existing set of nodes. VW points have still not been
        % inserted into the nodeset, which we'll do now. Similar to
        % U-line processing, we have to insert boundary elements as well.
%         disp('V-lines');
        nodeSearchFilter = false(allocSize, 1);
        nodeSearchFilter(nodeSizeTillPrevFace(i)+1 : nodeSize) = true;
        for j=1:3
            e1 = E_F(i, mod(j+1, 3)+1);
            e2 = E_F(i, mod(j+2, 3)+1);
            
            minVI = max(min(VI_E{e1}), min(VI_E{e2}));
            maxVI = min(max(VI_E{e1}), max(VI_E{e2}));
            
            % idx into UIonE{e1} and UIonE{e2}
            [~, shift1] = ismember(minVI, VI_E{e1});
            if isempty(shift1)
                shift1 = 0;
            end
            idx1 = shift1 + (0:maxVI-minVI);
            [~, shift2] = ismember(minVI, VI_E{e2});
            if isempty(shift2)
                shift2 = 0;
            end
            idx2 = shift2 + (0:maxVI-minVI);
            
            for k=1:numel(idx1)
                % We have to interpolate b/w the vertices at idx1(k)
                % position in PUI_E{e1} and idx2(k) in PUI_E{e2}.
                % UI_E{e1/2} and VonUI_E{e1/2} follow the same indexing
                k1 = idx1(k);
                k2 = idx2(k);
                
                ePair = [e1 e2];
                
                P12 = [PVI_E{e1}(k1, :); PVI_E{e2}(k2, :)];
                
                v = VI_E{e1}(k1);   % V is the same for both

                U12 = [UonPVI_E{e1}(k1); UonPVI_E{e2}(k2)];
                W12 = [WonPVI_E{e1}(k1); WonPVI_E{e2}(k2)];
                
                % point with smaller U is at alpha = 0
                if abs(diff(U12)) > abs(diff(W12))
                    [~, idxSort] = sort(U12);
                    sorter = 1;
                else
                    [~, idxSort] = sort(W12);
                    sorter = 3;
                end
                U12 = U12(idxSort);
                P12 = P12(idxSort, :);
                ePair = ePair(idxSort);
                W12 = W12(idxSort);

                if W12(1) > W12(2)
                    rangeW = integersHalfOpen(W12(2), W12(1));
                    alphaW = (rangeW - W12(2)) / (W12(1) - W12(2));
                    alphaW = 1 - alphaW;
                else
                    rangeW = integersHalfOpen(W12(1), W12(2));
                    alphaW = (rangeW - W12(1)) / diff(W12);
                end
                
                localNodeSearchFilter = nodeSearchFilter &...
                    nodeVI==v & nodeUI>=0;
                
                rangeU = nodeUI(localNodeSearchFilter);
                nui = numel(rangeU);
                
                nodeIdx_u = find(localNodeSearchFilter);
                
                % sort the searched list
                [rangeU, sortIdx] = sort(rangeU);
                nodeIdx_u = nodeIdx_u(sortIdx);
                w_ui = nodeW(nodeIdx_u);
              
                nwi = numel(alphaW);
                wiNodes = (1 - alphaW)*P12(1, :) + alphaW*P12(2, :);
                u_wi = (1 - alphaW)*U12(1) + alphaW*U12(2);
                PVWI_F{i} = [PVWI_F{i}; wiNodes];
                V_VWI_F{i} = [V_VWI_F{i}; repmat(v, nwi, 1)];
                W_VWI_F{i} = [W_VWI_F{i}; rangeW];
                UonPVWI_F{i} = [UonPVWI_F{i}; u_wi];

                nodeIdx_w = nodeSize + (1:nwi)';
                node(nodeIdx_w, :) = wiNodes;
                nodeVI(nodeIdx_w) = v;
                nodeWI(nodeIdx_w) = rangeW;
                nodeU(nodeIdx_w) = u_wi;
                nodeV(nodeIdx_w) = v;
                nodeW(nodeIdx_w) = rangeW;
                F_node(nodeIdx_w) = i;
                
                % add elements on the boundary
                if isFBoundary(i)
                    % find existing nodes on edges
                    e1NodeIdx = find(E_node==ePair(1) & nodeVI==v);
                    e2NodeIdx = find(E_node==ePair(2) & nodeVI==v);

                    % combine vars for ui and wi and sort by u
                    u = [double(rangeU); u_wi];
                    w = [w_ui; double(rangeW)];
                    nodeIdx = [nodeIdx_u; nodeIdx_w];
                    n = nui + nwi;
                    if sorter==1
                        [~, sortIdx] = sort(u);
                    else
                        [~, sortIdx] = sort(w);
                    end
                    nodeIdx = nodeIdx(sortIdx);
                    
                    nodeIsBoundary(nodeIdx_w) = true;
                    
                    
                    if n
                        % add the internal elems
                        elem(elemSize + (1:n-1), :) = [nodeIdx(1:n-1), nodeIdx(2:n)];
                        elemSize = elemSize + n-1;
                        
                        % add two elems if there is any internal node
                        elem(elemSize + (1:2), :) = [...
                            e1NodeIdx, nodeIdx(1);
                            nodeIdx(end), e2NodeIdx];
                        elemSize = elemSize + 2;
                    else
                        % otherwise, just add the elem connecting the two boundary
                        % nodes
                        elem(elemSize + 1, :) = [e1NodeIdx, e2NodeIdx];
                        elemSize = elemSize + 1;
                    end
                end
                
                nodeSize = nodeSize + nwi;
            end 
        end
        
        % Finally for W-lines. Now, we don't need to insert any new nodes.
%         disp('W-lines');
        nodeSearchFilter(nodeSizeTillPrevFace(i)+1 : nodeSize) = true;
        for j=1:3
            e1 = E_F(i, mod(j+1, 3)+1);
            e2 = E_F(i, mod(j+2, 3)+1);
            
            minWI = max(min(WI_E{e1}), min(WI_E{e2}));
            maxWI = min(max(WI_E{e1}), max(WI_E{e2}));
            
            % idx into UIonE{e1} and UIonE{e2}
            [~, shift1] = ismember(minWI, WI_E{e1});
            if isempty(shift1)
                shift1 = 0;
            end
            idx1 = shift1 + (0:maxWI-minWI);
            [~, shift2] = ismember(minWI, WI_E{e2});
            if isempty(shift2)
                shift2 = 0;
            end
            idx2 = shift2 + (0:maxWI-minWI);
            
            for k=1:numel(idx1)
                % We have to interpolate b/w the vertices at idx1(k)
                % position in PUI_E{e1} and idx2(k) in PUI_E{e2}.
                % UI_E{e1/2} and VonUI_E{e1/2} follow the same indexing
                k1 = idx1(k);
                k2 = idx2(k);
                
                ePair = [e1 e2];
                                
                w = WI_E{e1}(k1);   % W is the same for both
                
                U12 = [UonPWI_E{e1}(k1); UonPWI_E{e2}(k2)];
                V12 = [VonPWI_E{e1}(k1); VonPWI_E{e2}(k2)];
                
                % point with smaller U or V is at alpha = 0
                if abs(diff(U12)) > abs(diff(V12))
                    [~, idxSort] = sort(U12);
                    sorter = 1;
                else
                    [~, idxSort] = sort(V12);
                    sorter = 2;
                end
                
                ePair = ePair(idxSort);
                
                localNodeSearchFilter = nodeSearchFilter &...
                    nodeWI==w & nodeUI>=0;
                
                rangeU = nodeUI(localNodeSearchFilter);
                nodeIdx_u = find(localNodeSearchFilter);
                nui = numel(rangeU);
                
                % sort the searched lists
                [rangeU, sortIdx] = sort(rangeU);
                nodeIdx_u = nodeIdx_u(sortIdx);
                v_ui = nodeV(nodeIdx_u);
                
                localNodeSearchFilter = nodeSearchFilter &...
                    nodeWI==w & nodeVI>=0;
                
                nodeIdx_v = find(localNodeSearchFilter);
                nvi = numel(nodeIdx_v);
                
                % sort V-lists by U-values as well
                [u_vi, sortIdx] = sort(nodeU(nodeIdx_v));
                nodeIdx_v = nodeIdx_v(sortIdx);
                rangeV = nodeV(nodeIdx_v);
                

                % add elements on the boundary
                if isFBoundary(i)
                    % find existing nodes on edges
                    e1NodeIdx = find(E_node==ePair(1) & nodeWI==w);
                    e2NodeIdx = find(E_node==ePair(2) & nodeWI==w);

                    % combine vars for ui and vi and sort by u
                    u = [double(rangeU); u_vi];
                    v = [v_ui; double(rangeV)];
                    nodeIdx = [nodeIdx_u; nodeIdx_v];
                    n = nui + nvi;
                    if sorter==1
                        [~, sortIdx] = sort(u);
                    else
                        [~, sortIdx] = sort(v);
                    end
                    nodeIdx = nodeIdx(sortIdx);
                    
                    
                    if n
                        % add the internal elems
                        elem(elemSize + (1:n-1), :) = [nodeIdx(1:n-1), nodeIdx(2:n)];
                        elemSize = elemSize + n-1;
                        
                        % add two elems if there is any internal node
                        elem(elemSize + (1:2), :) = [...
                            e1NodeIdx, nodeIdx(1);
                            nodeIdx(end), e2NodeIdx];
                        elemSize = elemSize + 2;
                    % otherwise, just add the elem connecting the two boundary
                    % nodes
                    else
                        elem(elemSize + 1, :) = [e1NodeIdx, e2NodeIdx];
                        elemSize = elemSize + 1;
                    end
                end
            end 
        end
    end
    
    numBNode = nodeSize;

    
    % State of final graph so far:
    % node-> all nodes lying on faces (tris)
    % elem-> all elems on the boundary
    % includes extra nodes/elems which need to be collapsed (later)
    
    % for each tet, we have the nodes on the bounding faces
    % first, we iterate over each (unordered) pair of these bounding faces
    % Then, for each such pair, we iterate over each (unordered) pair of
    % tex-coordinates. We move between the two faces along each integral
    % isoline of this pair, and find the intersections with the isoplane of
    % the third tex-coordinate.
    for i=1:nT
%         disp(['Tet # ' num2str(i)]);
        
        f_choice = nchoosek(1:4, 2);
        nodeSizeTillPrevTet(i) = nodeSize;
        % First, UV-lines. This will find ALL the internal nodes and
        % also all the elems lying on UV-lines
%         disp('UV-lines');
        for j=1:size(f_choice, 1)
            f1 = F_T(i, f_choice(j, 1));
            f2 = F_T(i, f_choice(j, 2));
            
            idx1 = (1:numel(U_UVI_F{f1}))';
            
            for k1 = 1:numel(idx1)
                u = U_UVI_F{f1}(k1);
                v = V_UVI_F{f1}(k1);
                
                k2 = find(U_UVI_F{f2} == u & V_UVI_F{f2} == v);
                
                if isempty(k2)
                    continue;
                end
                
                fPair = [f1 f2]';
                
                P12 = [PUVI_F{f1}(k1, :); PUVI_F{f2}(k2, :)];
                
                W12 = [WonPUVI_F{f1}(k1); WonPUVI_F{f2}(k2)];
                
                [W12, idxSort] = sort(W12);
                P12 = P12(idxSort, :);
                fPair = fPair(idxSort);
                
                range = integersOpen(W12(1), W12(2));
                alpha = (range - W12(1)) / diff(W12);
                
                n = numel(alpha);
                
                allNodes = (1 - alpha)*P12(1, :) + alpha*P12(2, :);
                
                % find existing nodes on faces
                f1NodeIdx = find(F_node==fPair(1) & nodeUI==u & nodeVI==v);
                f2NodeIdx = find(F_node==fPair(2) & nodeUI==u & nodeVI==v);
                
                % add the internal nodes first
                nodeIdx = nodeSize + (1:n)';
                node(nodeIdx, :) = allNodes;
                nodeUI(nodeIdx) = u;
                nodeVI(nodeIdx) = v;
                nodeWI(nodeIdx) = range;
                nodeU(nodeIdx) = u;
                nodeV(nodeIdx) = v;
                nodeW(nodeIdx) = range;
                T_node(nodeIdx) = i;
                                           
                
                if n
                    % add the internal elems
                    elem(elemSize + (1:n-1), :) = [nodeIdx(1:n-1), nodeIdx(2:n)];
                    elemSize = elemSize + n-1;
                
                    % add two elems if there is any internal node
                    elem(elemSize + (1:2), :) = [...
                        f1NodeIdx, nodeIdx(1);...
                        nodeIdx(end), f2NodeIdx];
                    elemSize = elemSize + 2;
                else
                    % otherwise, just add the elem connecting the two nodes
                    % on the faces f1 and f2
                    elem(elemSize + 1, :) = [f1NodeIdx, f2NodeIdx];
                    elemSize = elemSize + 1;
                end
                
                nodeSize = nodeSize + n;
            end
        end
        
        nodeSearchFilter = false(allocSize, 1);
        nodeSearchFilter(nodeSizeTillPrevTet(i)+1 : nodeSize) = true;
        % Then, trace out the elems on VW-lines. Note that the nodes
        % are already there.
%         disp('VW-lines');
        for j=1:size(f_choice, 1)
            f1 = F_T(i, f_choice(j, 1));
            f2 = F_T(i, f_choice(j, 2));
            
            idx1 = (1:numel(V_VWI_F{f1}))';
            
            for k1 = 1:numel(idx1)
                v = V_VWI_F{f1}(k1);
                w = W_VWI_F{f1}(k1);
                
                k2 = find(V_VWI_F{f2} == v & W_VWI_F{f2} == w, 1);
                
                if isempty(k2)
                    continue;
                end
                
                U12 = [UonPVWI_F{f1}(k1); UonPVWI_F{f2}(k2)];
                
                fPair = [f1 f2]';
                
                [~, sortIdx] = sort(U12);
                fPair = fPair(sortIdx);
                
                localNodeSearchFilter = nodeSearchFilter &...
                    nodeVI==v & nodeWI==w & nodeUI>=0;
                
                range = nodeUI(localNodeSearchFilter);
                n = numel(range);
                
                nodeIdx = find(localNodeSearchFilter);
                
                %sort the searched lists
                [~, sortIdx] = sort(range);
                nodeIdx = nodeIdx(sortIdx);
                
                % find existing nodes on faces
                f1NodeIdx = find(F_node==fPair(1) & nodeVI==v & nodeWI==w);
                f2NodeIdx = find(F_node==fPair(2) & nodeVI==v & nodeWI==w);
                                           
                
                if n
                    % add the internal elems
                    elem(elemSize + (1:n-1), :) = [nodeIdx(1:n-1), nodeIdx(2:n)];
                    elemSize = elemSize + n-1;
                
                    % add two elems if there is any internal node
                    elem(elemSize + (1:2), :) = [...
                        f1NodeIdx, nodeIdx(1);...
                        nodeIdx(end), f2NodeIdx];
                    elemSize = elemSize + 2;
                else
                    % otherwise, just add the elem connecting the two nodes
                    % on the faces f1 and f2
                    elem(elemSize + 1, :) = [f1NodeIdx, f2NodeIdx];
                    elemSize = elemSize + 1;
                end                
            end
        end
        
        % Finally, trace out the elems on WU-lines. This works exactly the
        % same way as tracing VW-lines.
%         disp('WU-lines');
        for j=1:size(f_choice, 1)
            f1 = F_T(i, f_choice(j, 1));
            f2 = F_T(i, f_choice(j, 2));
            
            idx1 = (1:numel(W_WUI_F{f1}))';
            
            for k1 = 1:numel(idx1)
                w = W_WUI_F{f1}(k1);
                u = U_WUI_F{f1}(k1);
                
                k2 = find(W_WUI_F{f2} == w & U_WUI_F{f2} == u, 1);
                
                if isempty(k2)
                    continue;
                end
                
                V12 = [VonPWUI_F{f1}(k1); VonPWUI_F{f2}(k2)];
                
                fPair = [f1 f2]';
                
                [~, sortIdx] = sort(V12);
                fPair = fPair(sortIdx);
                
                localNodeSearchFilter = nodeSearchFilter &...
                    nodeWI==w & nodeUI==u & nodeVI>=0;
                
                range = nodeVI(localNodeSearchFilter);
                n = numel(range);
                
                nodeIdx = find(localNodeSearchFilter);
                
                %sort the searched lists
                [~, sortIdx] = sort(range);
                nodeIdx = nodeIdx(sortIdx);
                
                % find existing nodes on faces
                f1NodeIdx = find(F_node==fPair(1) & nodeWI==w & nodeUI==u);
                f2NodeIdx = find(F_node==fPair(2) & nodeWI==w & nodeUI==u);
                                           
                
                if n
                    % add the internal elems
                    elem(elemSize + (1:n-1), :) = [nodeIdx(1:n-1), nodeIdx(2:n)];
                    elemSize = elemSize + n-1;
                
                    % add two elems if there is any internal node
                    elem(elemSize + (1:2), :) = [...
                        f1NodeIdx, nodeIdx(1);...
                        nodeIdx(end), f2NodeIdx];
                    elemSize = elemSize + 2;
                else
                    % otherwise, just add the elem connecting the two nodes
                    % on the faces f1 and f2
                    elem(elemSize + 1, :) = [f1NodeIdx, f2NodeIdx];
                    elemSize = elemSize + 1;
                end                
            end
        end
    end
    
    
    node(nodeSize+1:end, :) = [];
    elem(elemSize+1:end, :) = [];
    nodeIsBoundary(nodeSize+1:end) = [];
    nodeUI(nodeSize+1:end) = [];
    nodeVI(nodeSize+1:end) = [];
    nodeWI(nodeSize+1:end) = [];
    E_node(nodeSize+1:end) = [];
    F_node(nodeSize+1:end) = [];
    T_node(nodeSize+1:end) = [];

    clear nodeU nodeV nodeW;
    
    disp(['Pre-collapse graph size: (' num2str(nodeSize) ', ' num2str(elemSize) ')']);

    % And finally, we need to collapse edges. This is done by
    % following isolines emanating from each boundary vertex, and
    % collapsing all invalid vertices along a line in one go.
    
    newIdx_node = zeros(nodeSize, 1);
    isEBCurved_node = false(nodeSize, 1);
    isEBCurved_node(E_node > 0) = isEBCurved(E_node(E_node > 0));
        
    node2 = zeros(nodeSize, 3);
    nodeUI2 = -int32(ones(nodeSize, 1));
    nodeVI2 = -int32(ones(nodeSize, 1));
    nodeWI2 = -int32(ones(nodeSize, 1));
    nodeB2 = false(nodeSize, 1);
    nodeC2 = false(nodeSize, 1);
    oldIdx_node2 = zeros(nodeSize, 1);

    
    elem2 = zeros(elemSize, 2);
    elemC2 = false(elemSize, 1);
    
    traced_node = false(numBNode, 1);
    tracedInternalUV = false(nodeSize, 1);
    tracedInternalVW = tracedInternalUV;
    tracedInternalWU = tracedInternalUV;
    
    adjMat_node = adjacency_matrix(elem);
    
    nodeSize = 0;
    elemSize = 0;
    
    disp('Collapsing process started.');
    disp('Finding boundary nodes...');
    
    
    if ~DONOTTRACEBOUNDARY
        disp('And tracing boundary curves...');
    end
    
    % Let's assume a max. of 10000 curves, with a max. of 1000 points each
    maxCC = 10000;
    cc = 0;
    maxPC = 1000;
    % For each curves, store the sequence of points, uvw-params, and 
    % binary variables declaring if the curve is a boundry curve (b),
    % a feature curve (f), and a closed curve (c).
    curves(maxCC) = struct('pts', zeros(maxPC, 1), ...
                           'u', -1, 'v', -1, 'w', -1, ...
                           'b', 0, 'f', 0, 'c', 0);
    
    
    % Trace UB curves: elems on the boundary with constant U-values
    boundaryNodes = find(...
        nodeIsBoundary &... % nodes on the boundary
        nodeUI >= 0 &... % with integer U-coordinate
        (nodeVI >= 0 | nodeWI >= 0 |... % and either integer W- or V-coordinate
        (E_node > 0 & isEBCurved_node))... % or lying on a curved edge
    );
    untracedBoundaryNodes = true(size(boundaryNodes));
    
    
%     disp('UB-curves');
    while any(untracedBoundaryNodes)
%         disp('Tracing boundary connected component...');
        firstNode = boundaryNodes(find(untracedBoundaryNodes, 1));
        
        u = nodeUI(firstNode);
        
        cc = cc + 1;
        curves(cc) = curves(maxCC);
        curves(cc).u = u;
        curves(cc).b = 1;
        pc = 0;
        
        untracedBoundaryNodes(boundaryNodes == firstNode) = false;
        
        if newIdx_node(firstNode)==0
            nodeSize = nodeSize + 1;
            node2(nodeSize, :) = node(firstNode, :);
            nodeB2(nodeSize) = true;
            if nodeUI(firstNode) >= 0
                nodeUI2(nodeSize) = nodeUI(firstNode);
            end
            if nodeVI(firstNode) >= 0
                nodeVI2(nodeSize) = nodeVI(firstNode);
            end
            if nodeWI(firstNode) >= 0
                nodeWI2(nodeSize) = nodeWI(firstNode);
            end
            
            firstNodeIdx = nodeSize;
            oldIdx_node2(nodeSize) = firstNode;
            newIdx_node(firstNode) = nodeSize;
        else
            firstNodeIdx = newIdx_node(firstNode);
        end
        
        curNode = firstNode;
        curNodeIdx = firstNodeIdx;
        prevNode = -1;
        
        pc = pc + 1;
        curves(cc).pts(pc) = firstNodeIdx;
        
        while (curNode ~= firstNode || prevNode==-1)
%             disp(curNode);
            adj = find(adjMat_node(:, curNode));
            nextNode = ...
                nodeUI(adj) == u &...   % follow the U-line
                nodeIsBoundary(adj) &...% stay on the boundary
                adj ~= prevNode &...        % don't go backwards, only fwd
                (F_node(adj) == F_node(curNode) |... % stay on the same face
                E_node(curNode) > 0 |... % or if currently on an edge
                E_node(adj)>0);   % or moving to an edge
            
            nextNode = adj(nextNode);
            
            if (prevNode == -1)
                nextNode = nextNode(1);
            end
            
            if numel(nextNode) ~= 1
                disp('Weird...');
            end
            
            if nodeVI(nextNode) >= 0 || nodeWI(nextNode) >= 0 || ...    % usual case: node lies on an integer isoline
                    (E_node(nextNode) > 0 && isEBCurved(E_node(nextNode)))  % to trace boundary geometry better, include
                                                                            % nodes on curved parts of the boundary
                untracedBoundaryNodes(boundaryNodes == nextNode) = false;
                if newIdx_node(nextNode) == 0
                    nodeSize = nodeSize + 1;
                    node2(nodeSize, :) = node(nextNode, :);
                    nodeB2(nodeSize) = true;
                    if nodeUI(nextNode) >= 0
                        nodeUI2(nodeSize) = nodeUI(nextNode);
                    end
                    if nodeVI(nextNode) >= 0
                        nodeVI2(nodeSize) = nodeVI(nextNode);
                    end
                    if nodeWI(nextNode) >= 0
                        nodeWI2(nodeSize) = nodeWI(nextNode);
                    end
                    oldIdx_node2(nodeSize) = nextNode;
                    nextNodeIdx = nodeSize;
                    newIdx_node(nextNode) = nodeSize;
                else
                    nextNodeIdx = newIdx_node(nextNode);
                end
                
                if ~DONOTTRACEBOUNDARY
                    elem2(elemSize + 1, :) = [curNodeIdx, nextNodeIdx];
                    pc = pc + 1;
                    curves(cc).pts(pc) = nextNodeIdx;
					if E_node(nextNode) > 0 && isEBCurved(E_node(nextNode))
						curves(cc).pts(pc+1:end) = [];
                        cc = cc + 1;
						curves(cc) = curves(maxCC);
						curves(cc).u = u;
						curves(cc).b = 1;
						pc = 1;
						curves(cc).pts(pc) = nextNodeIdx;
					end
                    elemSize = elemSize + 1;
                end
                curNodeIdx = nextNodeIdx;
            end
            prevNode = curNode;
            curNode = nextNode;
        end
        
        curves(cc).pts(pc+1:end) = [];
    end
    
    boundaryNodes = find(...
        nodeIsBoundary &... % nodes on the boundary
        nodeVI >= 0 &... % with integer V-coordinate
        (nodeUI >= 0 | nodeWI >= 0 |... % and either integer U- or W-coordinate
        (E_node > 0 & isEBCurved_node))... % or lying on a curved edge)
    );
    untracedBoundaryNodes = true(size(boundaryNodes));
%     disp('VB-curves');
    while any(untracedBoundaryNodes)
%         disp('Tracing boundary connected component...');
        
        firstNode = boundaryNodes(find(untracedBoundaryNodes, 1));
        
        v = nodeVI(firstNode);
        
        cc = cc + 1;
        curves(cc) = curves(maxCC);
        curves(cc).v = v;
        curves(cc).b = 1;
        pc = 0;
        
        untracedBoundaryNodes(boundaryNodes == firstNode) = false;
        
        if newIdx_node(firstNode)==0
            nodeSize = nodeSize + 1;
            node2(nodeSize, :) = node(firstNode, :);
            nodeB2(nodeSize) = true;
            if nodeUI(firstNode) >= 0
                nodeUI2(nodeSize) = nodeUI(firstNode);
            end
            if nodeVI(firstNode) >= 0
                nodeVI2(nodeSize) = nodeVI(firstNode);
            end
            if nodeWI(firstNode) >= 0
                nodeWI2(nodeSize) = nodeWI(firstNode);
            end
            firstNodeIdx = nodeSize;
            oldIdx_node2(nodeSize) = firstNode;
            newIdx_node(firstNode) = nodeSize;
        else
            firstNodeIdx = newIdx_node(firstNode);
        end

        curNode = firstNode;
        curNodeIdx = firstNodeIdx;
        prevNode = -1;
        
        pc = pc + 1;
        curves(cc).pts(pc) = firstNodeIdx;
        
        while (curNode ~= firstNode || prevNode==-1)
%             disp(curNode);
            adj = find(adjMat_node(:, curNode));
            nextNode = ...
                nodeVI(adj) == v &...   % follow the U-line
                nodeIsBoundary(adj) &...% stay on the boundary
                adj ~= prevNode &...        % don't go backwards, only fwd
                (F_node(adj) == F_node(curNode) |... % stay on the same face
                E_node(curNode) > 0 |... % or if currently on an edge
                E_node(adj)>0);   % or moving to an edge
            
            nextNode = adj(nextNode);
            
            if prevNode == -1
                nextNode = nextNode(1);
            end
            
            if numel(nextNode) ~= 1
                disp('Weird...');
            end
            
            if nodeUI(nextNode) >= 0 || nodeWI(nextNode) >= 0 || ...    % usual case: node lies on an integer isoline
                    (E_node(nextNode) > 0 && isEBCurved(E_node(nextNode)))  % to trace boundary geometry better, include
                                                                            % nodes on curved parts of the boundary
                untracedBoundaryNodes(boundaryNodes == nextNode) = false;
                if newIdx_node(nextNode) == 0
                    nodeSize = nodeSize + 1;
                    node2(nodeSize, :) = node(nextNode, :);
                    nodeB2(nodeSize) = true;
                    if nodeUI(nextNode) >= 0
                        nodeUI2(nodeSize) = nodeUI(nextNode);
                    end
                    if nodeVI(nextNode) >= 0
                        nodeVI2(nodeSize) = nodeVI(nextNode);
                    end
                    if nodeWI(nextNode) >= 0
                        nodeWI2(nodeSize) = nodeWI(nextNode);
                    end
                    oldIdx_node2(nodeSize) = nextNode;
                    nextNodeIdx = nodeSize;
                    newIdx_node(nextNode) = nodeSize;
                else
                    nextNodeIdx = newIdx_node(nextNode);
                end
                
                if ~DONOTTRACEBOUNDARY
                    elem2(elemSize + 1, :) = [curNodeIdx, nextNodeIdx];
                    pc = pc + 1;
                    curves(cc).pts(pc) = nextNodeIdx;
					if E_node(nextNode) > 0 && isEBCurved(E_node(nextNode))
						curves(cc).pts(pc+1:end) = [];
                        cc = cc + 1;
						curves(cc) = curves(maxCC);
						curves(cc).v = v;
						curves(cc).b = 1;
						pc = 1;
						curves(cc).pts(pc) = nextNodeIdx;
					end
                    elemSize = elemSize + 1;
                end
                curNodeIdx = nextNodeIdx;
            end
            prevNode = curNode;
            curNode = nextNode;
        end
        
        curves(cc).pts(pc+1:end) = [];
    end

    boundaryNodes = find(...
        nodeIsBoundary &... % nodes on the boundary
        nodeWI >= 0 &... % with integer W-coordinate
        (nodeUI >= 0 | nodeVI >= 0 |... % and either integer U- or V-coordinate
        (E_node > 0 & isEBCurved_node))... % or lying on a curved edge
    );
    untracedBoundaryNodes = true(size(boundaryNodes));
%     disp('WB-curves');
    while any(untracedBoundaryNodes)
%         disp('Tracing boundary connected component...');
        
        firstNode = boundaryNodes(find(untracedBoundaryNodes, 1));
        
        w = nodeWI(firstNode);
        
        cc = cc + 1;
        curves(cc) = curves(maxCC);
        curves(cc).w = w;
        curves(cc).b = 1;
        pc = 0;
        
        untracedBoundaryNodes(boundaryNodes == firstNode) = false;
        
        if newIdx_node(firstNode)==0
            nodeSize = nodeSize + 1;
            node2(nodeSize, :) = node(firstNode, :);
            nodeB2(nodeSize) = true;
            if nodeUI(firstNode) >= 0
                nodeUI2(nodeSize) = nodeUI(firstNode);
            end
            if nodeVI(firstNode) >= 0
                nodeVI2(nodeSize) = nodeVI(firstNode);
            end
            if nodeWI(firstNode) >= 0
                nodeWI2(nodeSize) = nodeWI(firstNode);
            end
            firstNodeIdx = nodeSize;
            oldIdx_node2(nodeSize) = firstNode;
            newIdx_node(firstNode) = nodeSize;
        else
            firstNodeIdx = newIdx_node(firstNode);
        end
        
        curNode = firstNode;
        curNodeIdx = firstNodeIdx;
        prevNode = -1;
        
        pc = pc + 1;
        curves(cc).pts(pc) = firstNodeIdx;
        
        while (curNode ~= firstNode || prevNode==-1)
%             disp([curNode curNodeIdx nodeUI(curNode) nodeVI(curNode) E_node(curNode)]);
            adj = find(adjMat_node(:, curNode));
            nextNode = ...
                nodeWI(adj) == w &...   % follow the W-line
                nodeIsBoundary(adj) &...% stay on the boundary
                adj ~= prevNode &...        % don't go backwards, only fwd
                (F_node(adj) == F_node(curNode) |... % stay on the same face
                E_node(curNode) > 0 |... % or if currently on an edge
                E_node(adj)>0);   % or moving to an edge
            
            nextNode = adj(nextNode);
            
            if prevNode == -1
                nextNode = nextNode(1);
            end
            
            if numel(nextNode) ~= 1
                disp('Weird...');
            end
            
            if nodeUI(nextNode) >= 0 || nodeVI(nextNode) >= 0 || ...    % usual case: node lies on an integer isoline
                    (E_node(nextNode) > 0 && isEBCurved(E_node(nextNode)))  % to trace boundary geometry better, include
                                                                            % nodes on curved parts of the boundary
                untracedBoundaryNodes(boundaryNodes == nextNode) = false;
                if newIdx_node(nextNode) == 0
                    nodeSize = nodeSize + 1;
                    node2(nodeSize, :) = node(nextNode, :);
                    nodeB2(nodeSize) = true;
                    if nodeUI(nextNode) >= 0
                        nodeUI2(nodeSize) = nodeUI(nextNode);
                    end
                    if nodeVI(nextNode) >= 0
                        nodeVI2(nodeSize) = nodeVI(nextNode);
                    end
                    if nodeWI(nextNode) >= 0
                        nodeWI2(nodeSize) = nodeWI(nextNode);
                    end
                    oldIdx_node2(nodeSize) = nextNode;
                    nextNodeIdx = nodeSize;
                    newIdx_node(nextNode) = nodeSize;
                else
                    nextNodeIdx = newIdx_node(nextNode);
                end
                
                if ~DONOTTRACEBOUNDARY
                    elem2(elemSize + 1, :) = [curNodeIdx, nextNodeIdx];
                    pc = pc + 1;
                    curves(cc).pts(pc) = nextNodeIdx;
					if E_node(nextNode) > 0 && isEBCurved(E_node(nextNode))
                        curves(cc).pts(pc+1:end) = [];
						cc = cc + 1;
						curves(cc) = curves(maxCC);
						curves(cc).w = w;
						curves(cc).b = 1;
						pc = 1;
						curves(cc).pts(pc) = nextNodeIdx;
					end
                    elemSize = elemSize + 1;
                end
                curNodeIdx = nextNodeIdx;
            end
            prevNode = curNode;
            curNode = nextNode;
        end
        
        curves(cc).pts(pc+1:end) = [];
    end
    
    numBoundaryNode = nodeSize;
    
    disp(['# of boundary nodes: ', num2str(numBoundaryNode)]);

    if ~DONOTTRACEINTERNAL
        disp('Tracing internal curves now...');
        % Now collapse the internal edges
        for i=1:size(node, 1)
            if ~nodeIsBoundary(i) || ...    % lies on the domain boundary
                traced_node(i) || ...       % not traced already
                ... % lies on intersection of 2 integral isocurves (3 sub-conditions)
                (nodeUI(i) <= 0 && nodeVI(i) <= 0) ||... % U or V integer
                (nodeVI(i) <= 0 && nodeWI(i) <= 0) ||... % V or W integer
                (nodeWI(i) <= 0 && nodeUI(i) <= 0)       % W or U integer
                continue;
            end

            u = nodeUI(i);
            v = nodeVI(i);
            w = nodeWI(i);
            
            cc = cc + 1;
            curves(cc) = curves(maxCC);
            curves(cc).u = u;
            curves(cc).v = v;
            curves(cc).w = w;
                
            pc = 0;

            traced_node(i) = true;

            prevNode = 0;
            curNode = i;
            curNodeIdx = newIdx_node(curNode);
            
            pc = pc + 1;
            curves(cc).pts(pc) = curNodeIdx;

    %         disp(curNode);

            if u >= 0 && v >= 0
                while curNode == i || ~nodeIsBoundary(curNode)
                    adj = find(adjMat_node(:, curNode));
                    nextNode = nodeUI(adj)==u & nodeVI(adj)==v &...    % follow the line
                        adj ~= prevNode;   % don't go back, only forward

                    nextNode = adj(nextNode);

                    if numel(nextNode) ~= 1
                        disp('Weird...');
                    end

                    % Only points on W-isocountour or boundary should be 
                    % added to the final graph
                    if nodeWI(nextNode)>=0 || nodeIsBoundary(nextNode)
                        nextNodeIdx = newIdx_node(nextNode);
                        if nextNodeIdx == 0
                            nodeSize = nodeSize + 1;
                            node2(nodeSize, :) = node(nextNode, :);
                            nodeUI2(nodeSize) = nodeUI(nextNode);
                            nodeVI2(nodeSize) = nodeVI(nextNode);
                            if nodeWI(nextNode) >= 0
                                nodeWI2(nodeSize) = nodeWI(nextNode);
                            end
                            oldIdx_node2(nodeSize) = nextNode;
                            nextNodeIdx = nodeSize;
                            newIdx_node(nextNode) = nodeSize;
                        end
%                         if ~DONOTTRACEINTERNAL
                            elem2(elemSize + 1, :) = [curNodeIdx, nextNodeIdx];
                            tracedInternalUV([oldIdx_node2(curNodeIdx), nextNode]) = true;
                            pc = pc + 1;
                            curves(cc).pts(pc) = nextNodeIdx;
                            elemSize = elemSize + 1;
%                         end
                        curNodeIdx = nextNodeIdx;
                    end
                    prevNode = curNode;
                    curNode = nextNode;
                end
            elseif v >= 0 && w >= 0
                while curNode == i || ~nodeIsBoundary(curNode)
                    adj = find(adjMat_node(:, curNode));
                    nextNode = nodeVI(adj)==v & nodeWI(adj)==w &...    % follow the line
                        adj ~= prevNode;   % don't go back, only forward

                    nextNode = adj(nextNode);

                    if numel(nextNode) ~= 1
                        disp('Weird...');
                    end

                    % Only UVW-intersections and boundary points should be 
                    % added to the final graph
                    if nodeUI(nextNode)>=0 || nodeIsBoundary(nextNode)
                        nextNodeIdx = newIdx_node(nextNode);
                        if nextNodeIdx == 0
                            nodeSize = nodeSize + 1;
                            node2(nodeSize, :) = node(nextNode, :);
                            nodeVI2(nodeSize) = nodeVI(nextNode);
                            nodeWI2(nodeSize) = nodeWI(nextNode);
                            if nodeUI(nextNode) >= 0
                                nodeUI2(nodeSize) = nodeUI(nextNode);
                            end
                            oldIdx_node2(nodeSize) = nextNode;
                            nextNodeIdx = nodeSize;
                            newIdx_node(nextNode) = nodeSize;
                        end
%                         if ~DONOTTRACEINTERNAL
                            elem2(elemSize + 1, :) = [curNodeIdx, nextNodeIdx];
                            tracedInternalVW([oldIdx_node2(curNodeIdx), nextNode]) = true;
                            pc = pc + 1;
                            curves(cc).pts(pc) = nextNodeIdx;
                            elemSize = elemSize + 1;
%                         end
                        curNodeIdx = nextNodeIdx;
                    end
                    prevNode = curNode;
                    curNode = nextNode;
                end
            elseif w >= 0 && u >= 0
                while curNode == i || ~nodeIsBoundary(curNode)
                    adj = find(adjMat_node(:, curNode));
                    nextNode = nodeWI(adj)==w & nodeUI(adj)==u &...    % follow the line
                        adj ~= prevNode;   % don't go back, only forward

                    nextNode = adj(nextNode);

                    if numel(nextNode) ~= 1
                        disp('Weird...');
                    end

                    % Only UVW-intersections and boundary points should be 
                    % added to the final graph
                    if nodeVI(nextNode)>=0 || nodeIsBoundary(nextNode)
                        nextNodeIdx = newIdx_node(nextNode);
                        if nextNodeIdx == 0
                            nodeSize = nodeSize + 1;
                            node2(nodeSize, :) = node(nextNode, :);
                            nodeWI2(nodeSize) = nodeWI(nextNode);
                            nodeUI2(nodeSize) = nodeUI(nextNode);
                            if nodeVI(nextNode) >= 0
                                nodeVI2(nodeSize) = nodeVI(nextNode);
                            end
                            oldIdx_node2(nodeSize) = nextNode;
                            nextNodeIdx = nodeSize;
                            newIdx_node(nextNode) = nodeSize;
                        end
%                         if ~DONOTTRACEINTERNAL
                            elem2(elemSize + 1, :) = [curNodeIdx, nextNodeIdx];
                            tracedInternalWU([oldIdx_node2(curNodeIdx), nextNode]) = true;
                            pc = pc + 1;
                            curves(cc).pts(pc) = nextNodeIdx;
                            elemSize = elemSize + 1;
%                         end
                        curNodeIdx = nextNodeIdx;
                    end
                    prevNode = curNode;
                    curNode = nextNode;
                end
            end

            traced_node(curNode) = true;
            
            curves(cc).pts(pc+1:end) = [];
        end
        
        % Trace closed isocurves now: The method will be fairly similar to
        % tracing curves on the boundary.
        deg = accumarray(reshape(elem2(1:elemSize, :), [2*elemSize, 1]), 1);
        deg = [deg; zeros(nodeSize - numel(deg), 1)];
        
        isNodeUntraced = nodeUI2(1:nodeSize)>=0 & nodeVI2(1:nodeSize)>=0 & nodeWI2(1:nodeSize)>=0 & deg < 6;
        untracedNodeIdx = find(isNodeUntraced);
        numUntraced = numel(untracedNodeIdx);
        untracedNodeIdx_old = oldIdx_node2(untracedNodeIdx);
        uvUntraced = 2*int32(~tracedInternalUV(untracedNodeIdx_old(1:numUntraced)));
        vwUntraced = 2*int32(~tracedInternalVW(untracedNodeIdx_old(1:numUntraced)));
        wuUntraced = 2*int32(~tracedInternalWU(untracedNodeIdx_old(1:numUntraced)));
        

        i = 1;
        while(i <= numUntraced)
            if ~uvUntraced(i) && ~vwUntraced(i) && ~wuUntraced(i)
                i = i+1;
                continue
            end
            
            firstNode = untracedNodeIdx_old(i);
            
            u = nodeUI(firstNode);
            v = nodeVI(firstNode);
            w = nodeWI(firstNode);
            
            cc = cc + 1;
            curves(cc) = curves(maxCC);
            curves(cc).u = u;
            curves(cc).v = v;
            curves(cc).w = w;
            curves(cc).c = 1;
            
            pc = 1;
            curves(cc).pts(pc) = curNodeIdx;
            
            prevNode = -1;
            curNode = firstNode;
            curNodeIdx = newIdx_node(curNode);

            if uvUntraced(i)
%                 disp(['Tracing UV for node#', num2str(i)]);
%                 disp(['uvUntraced(i) = ', num2str(uvUntraced(i))]);
                while curNode ~= firstNode || prevNode == -1
                    adj = find(adjMat_node(:, curNode));
                    nextNode = nodeUI(adj)==u & nodeVI(adj)==v &...    % follow the line
                        adj ~= prevNode;   % don't go back, only forward

                    nextNode = adj(nextNode);
                    
                    if prevNode == -1
                        nextNode = nextNode(1);
                    end
                    
                    if numel(nextNode) == 0
                        disp('Curve ended!');
                    end

                    if numel(nextNode) ~= 1
                        disp('Weird...');
                    end

                    % Only UVW-intersections should be 
                    % added to the final graph
                    if nodeWI(nextNode)>=0
                        nextNodeIdx = newIdx_node(nextNode);
                        elem2(elemSize + 1, :) = [curNodeIdx, nextNodeIdx];
                        uvUntraced(untracedNodeIdx == curNodeIdx) = uvUntraced(untracedNodeIdx == curNodeIdx) - 1;
                        uvUntraced(untracedNodeIdx == nextNodeIdx) = uvUntraced(untracedNodeIdx == nextNodeIdx) - 1;
                        pc = pc + 1;
                        curves(cc).pts(pc) = nextNodeIdx;
                        elemSize = elemSize + 1;
                        curNodeIdx = nextNodeIdx;
                    end
                    prevNode = curNode;
                    curNode = nextNode;
                end
            elseif vwUntraced(i)
%                 disp(['Tracing VW for node#', num2str(i)]);
%                 disp(['vwUntraced(i) = ', num2str(vwUntraced(i))]);
                while curNode ~= firstNode || prevNode == -1
                    adj = find(adjMat_node(:, curNode));
                    nextNode = nodeVI(adj)==v & nodeWI(adj)==w &...    % follow the line
                        adj ~= prevNode;   % don't go back, only forward

                    nextNode = adj(nextNode);
                    
                    if prevNode == -1
                        nextNode = nextNode(1);
                    end

                    if numel(nextNode) == 0
                        disp('Curve ended!');
                    end
                    
                    if numel(nextNode) ~= 1
                        disp('Weird...');
                    end

                    % Only UVW-intersections should be 
                    % added to the final graph
                    if nodeUI(nextNode)>=0
                        nextNodeIdx = newIdx_node(nextNode);
                        elem2(elemSize + 1, :) = [curNodeIdx, nextNodeIdx];
                        vwUntraced(untracedNodeIdx == curNodeIdx) = vwUntraced(untracedNodeIdx == curNodeIdx) - 1;
                        vwUntraced(untracedNodeIdx == nextNodeIdx) = vwUntraced(untracedNodeIdx == nextNodeIdx) - 1;
                        pc = pc + 1;
                        curves(cc).pts(pc) = nextNodeIdx;
                        elemSize = elemSize + 1;
                        curNodeIdx = nextNodeIdx;
                    end
                    prevNode = curNode;
                    curNode = nextNode;
                end
            elseif wuUntraced(i)
%                 disp(['Tracing WU for node#', num2str(i)]);
%                 disp(['wuUntraced(i) = ', num2str(wuUntraced(i))]);
                while curNode ~= firstNode || prevNode == -1
                    adj = find(adjMat_node(:, curNode));
                    nextNode = nodeWI(adj)==w & nodeUI(adj)==u &...    % follow the line
                        adj ~= prevNode;   % don't go back, only forward

                    nextNode = adj(nextNode);

                    if prevNode == -1
                        nextNode = nextNode(1);
                    end

                    if numel(nextNode) == 0
                        disp('Curve ended!');
                    end
                    
                    if numel(nextNode) ~= 1
                        disp('Weird...');
                    end

                    % Only UVW-intersections should be 
                    % added to the final graph
                    if nodeVI(nextNode)>=0
                        nextNodeIdx = newIdx_node(nextNode);
                        elem2(elemSize + 1, :) = [curNodeIdx, nextNodeIdx];
                        wuUntraced(untracedNodeIdx == curNodeIdx) = wuUntraced(untracedNodeIdx == curNodeIdx) - 1;
                        wuUntraced(untracedNodeIdx == nextNodeIdx) = wuUntraced(untracedNodeIdx == nextNodeIdx) - 1;
                        pc = pc + 1;
                        curves(cc).pts(pc) = nextNodeIdx;
                        elemSize = elemSize + 1;
                        curNodeIdx = nextNodeIdx;
                    end
                    prevNode = curNode;
                    curNode = nextNode;
                end
            end
            
            curves(cc).pts(pc+1:end) = [];
        end
    end
    
    curves(cc+1:end) = [];
    
    % trace out edge features
    elemBC = zeros(100000, 2);
    nodeBC = zeros(50000, 3);
    nodeBCSize = 0;
    elemBCSize = 0;
    

    ccBC = 0;

    % For each curves, store the sequence of points, uvw-params, and 
    % binary variables declaring if the curve is a boundry curve (b),
    % a feature curve (f), and a closed curve (c).
    curvesBC(maxCC) = struct('pts', zeros(maxPC, 1), ...
                           'u', -1, 'v', -1, 'w', -1, ...
                           'b', 1, 'f', 1, 'c', 0);

    if ~DONOTTRACEBOUNDARY
        % find degree of vertices in the subgraph (P, EBCurved)
        EBCurved = E(isEBCurved, :);
        D = accumarray(EBCurved(:), 1);
        D = [D; zeros(nP - numel(D), 1)];
        
        % while any untraced feature curves are left
        while any(D)
            
            ccBC = ccBC + 1;
            curvesBC(ccBC) = curvesBC(maxCC);
            pc = 0;
            
            % get first point with D(p) = 1, which indicates that it is at
            % the end of a feature curve
            p1Idx = find(D==1, 1);
            % maybe we only have closed feature curves left
            if isempty(p1Idx)
                p1Idx = find(D>0 & isPBCurved, 1);
            end
            if isempty(p1Idx)
                p1Idx = find(D>0, 1);
            end
            p1 = P(p1Idx, :);
            
            nodeBCSize = nodeBCSize + 1;
            nodeBC(nodeBCSize, :) = p1;
            firstNodeIdx = nodeSize + nodeBCSize;
            
            curNodeIdx = firstNodeIdx;
            
            pc = pc + 1;
            curvesBC(ccBC).pts(pc) = firstNodeIdx;
            
            % find a curved edge incident on p
            eIdx = find(isEBCurved & (E(:, 1) == p1Idx | E(:, 2) == p1Idx), 1);

            while ~isempty(eIdx)
                % update all data structures
                isEBCurved(eIdx) = false;
                e = E(eIdx, :);
                p2Idx = e(find(e~=p1Idx, 1));
                D(p1Idx) = D(p1Idx) - 1;
                D(p2Idx) = D(p2Idx) - 1;
                p2 = P(p2Idx, :);
                
                nodeIdx_e = newIdx_node(E_node == eIdx);
                if ~isempty(nodeIdx_e)
                    node_e = node2(nodeIdx_e, :);
                    % get extants of these nodes
                    extants = max(node_e) - min(node_e);
                    [~, sortAxis] = max(extants);
                    [~, sortIdx] = sort(node_e(:, sortAxis));
                    if p1(sortAxis) < p2(sortAxis)
                        nodeIdx_e = nodeIdx_e(sortIdx);
                    else
                        nodeIdx_e = nodeIdx_e(flip(sortIdx));
                    end
                    
                    % get the elements lying on this edge, and the one
                    % connecting curNodeIdx to the first node on this edge
                    elemBC(elemBCSize + (1:numel(nodeIdx_e)), :) = ...
                        [[curNodeIdx; nodeIdx_e(1:end-1)], nodeIdx_e];
                    curvesBC(ccBC).pts(pc + (1:numel(nodeIdx_e))) = nodeIdx_e;
                    pc = pc + numel(nodeIdx_e);
                    elemBCSize = elemBCSize + numel(nodeIdx_e);
                    
                    curNodeIdx = nodeIdx_e(end);
                end
                
                % keep adding curved vertices found on the way
                % this can add the same vertex many times, but duplicate
                % nodes are easy to deal with here
                if isPBCurved(p2Idx)
                    nodeBCSize = nodeBCSize + 1;
                    nodeBC(nodeBCSize, :) = p2;
                    elemBCSize = elemBCSize + 1;
                    elemBC(elemBCSize, :) = ...
                        [curNodeIdx, nodeSize + nodeBCSize];
                    
                    pc = pc + 1;
                    curvesBC(ccBC).pts(pc) = nodeSize + nodeBCSize;
                    
                    % while we make sure no edges are repeated and thus the
                    % path thus not loop over itself, it may still touch
                    % the initial point since we do not ensure that
                    % vertices are not repeated. Therefore, end the current
                    % curve and start a new one if the vertex in
                    % consideration still has a curved edge incident on it
                    if D(p2Idx)
                        curvesBC(ccBC).pts(pc+1:end) = [];
                        pc = 1;
                        ccBC = ccBC + 1;
                        curvesBC(ccBC) = curvesBC(maxCC);
                        curvesBC(ccBC).pts(pc) = nodeSize + nodeBCSize;
                    end
                    curNodeIdx = nodeSize + nodeBCSize;
                end
                
                % update p1Idx and eIdx
                p1Idx = p2Idx;
                p1 = p2;
                eIdx = find(isEBCurved & (E(:, 1) == p1Idx | E(:, 2) == p1Idx), 1);
            end
            
            % then, connect curNodeIdx to the last vertex of this
            % feature curve, if ot done already
            if ~isPBCurved(p1Idx)
                nodeBCSize = nodeBCSize + 1;
                nodeBC(nodeBCSize, :) = p1;
                lastNodeIdx = nodeSize + nodeBCSize;
                elemBCSize = elemBCSize + 1;
                elemBC(elemBCSize, :) = [curNodeIdx, lastNodeIdx];
                
                pc = pc + 1;
                curvesBC(ccBC).pts(pc) = lastNodeIdx;
            end
            
            curvesBC(ccBC).pts(pc+1:end) = [];
        end
        
        nodeBC(nodeBCSize+1:end, :) = [];
        elemBC(elemBCSize+1:end, :) = [];
        curvesBC(ccBC+1:end) = [];
        
        % Finally, insert these newly traced boundary nodes and elements in
        % node2 and elem2. We need to ensure that all boundary nodes form
        % the top #NB rows of node2. The elements can be inserted anywhere.
        
        % first, insert at the end
        node2 = [node2(1:nodeSize, :); nodeBC];
        nodeB2 = [nodeB2(1:nodeSize); true(nodeBCSize, 1)];
        nodeC2 = [nodeC2(1:nodeSize); true(nodeBCSize, 1)];
        elem2 = [elem2(1:elemSize, :); elemBC];
        elemC2 = [elemC2(1:elemSize); true(elemBCSize, 1)];
        
        % now, move these new nodes to the top and modify elem2 accordingly
        idxInv = [ nodeBCSize+(1:nodeSize), (1:nodeBCSize) ]';
        idx = [nodeSize+(1:nodeBCSize), 1:nodeSize]';
        
        node2 = node2(idx, :);
        nodeB2 = nodeB2(idx);
        nodeC2 = nodeC2(idx);
        nodeUI2 = nodeUI2(idx);
        nodeVI2 = nodeVI2(idx);
        nodeWI2 = nodeWI2(idx);
        nodeSize = nodeSize + nodeBCSize;
        
        elem2 = idxInv(elem2);
        elemSize = elemSize + elemBCSize;
        
        curves = [curves, curvesBC];
        
        for cc=1:numel(curves)
            curves(cc).pts = idxInv(curves(cc).pts);
        end
%         numBoundaryNode = numBoundaryNode + nodeBCSize;
        
%         disp(['# of boundary nodes (possibly) increased to ', ...
%             num2str(numBoundaryNode)]);
    end
    
    elem2(elemSize+1 : end, :) = [];
    elemC2(elemSize+1:end) = [];
    
    node2(nodeSize+1 : end, :) = [];
    nodeB2(nodeSize+1:end) = [];
    nodeC2(nodeSize+1:end) = [];
    nodeUI2(nodeSize+1:end) = [];
    nodeVI2(nodeSize+1:end) = [];
    nodeWI2(nodeSize+1:end) = [];
    
    % collapse duplicate nodes
    [node2, SVI, SVJ] = remove_duplicate_vertices(node2);
    nodeB2 = nodeB2(SVI);
    nodeC2 = nodeC2(SVI);
    nodeUI2 = nodeUI2(SVI);
    nodeVI2 = nodeVI2(SVI);
    nodeWI2 = nodeWI2(SVI);
    elem2 = SVJ(elem2);
    elemC2(elem2(:, 1)==elem2(:, 2)) = [];
    elem2(elem2(:, 1)==elem2(:, 2), :) = [];
    
    for cc=1:numel(curves)
        curves(cc).pts = SVJ(curves(cc).pts);
        if curves(cc).pts(1) == curves(cc).pts(end)
            curves(cc).c = 1;
        else
            curves(cc).c = 0;
        end
    end
    
	
    % stupid hack to remove disconnected elements
    d = accumarray(elem2(:), 1);
    elemC2(d(elem2(:, 1))==1) = [];
    elem2(d(elem2(:, 1))==1, :) = [];
    elemC2(d(elem2(:, 2))==1) = [];
    elem2(d(elem2(:, 2))==1, :) = [];
    elemSize = size(elem2, 1);
    
    for cc=1:numel(curves)
        if d(curves(cc).pts(1))==1
            curves(cc).pts(1) = [];
        end
        if d(curves(cc).pts(end))==1
            curves(cc).pts(end) = [];
        end
    end
    
    node = node2;
    elem = elem2;
    
    disp(['Post-collapse graph size: (' num2str(nodeSize) ', ' num2str(elemSize) ')']);
    
    dataOut = dataIn;
    dataOut.Node = node;
    dataOut.Elem = elem;
    dataOut.NodeB = nodeB2;
    dataOut.NodeC = nodeC2;
    dataOut.NodeU = nodeUI2;
    dataOut.NodeV = nodeVI2;
    dataOut.NodeW = nodeWI2;
    dataOut.ElemC = elemC2;
    dataOut.Curves = curves;
end

function r = integersOpen(a, b, isBoundary)
    r = (ceil(a):floor(b))';
end

function r = integersHalfOpen(a, b, isBoundary)
    r = ((floor(a) + 1) : floor(b))';
end