function dataOut = fitFramesToData3D(V, T, stress)

    tol = sqrt(eps);
    
    nT = size(T,1);
    nV = size(V, 1);
    w = zeros(size(V,1), 3); %axis angle rep 
    L = cotmatrix(V,T); %using good old laplacian as a smoothing energy (kind of like static elasticity)

    H = [L'*L 0*L 0*L;0*L L'*L 0*L; 0*L 0*L L'*L];

    %mesh indexing sparse matrix for gathers
    FI = repmat(1:nT, size(T,2),1);
    FT = T';
    V2F = sparse(FI(:), FT(:), ones(numel(T),1), nT, nV);
	
    %setup stress tensors
    %0 xx 1 yy 2 zz 3 xy 4 yz 5 xz
    sxx = stress(1:nT);
    syy = stress(nT + (1:nT));
    szz = stress(2*nT + (1:nT));
    sxy = stress(3*nT + (1:nT));
    syz = stress(4*nT + (1:nT));
    sxz = stress(5*nT + (1:nT));

    S1 = zeros(nT, 3);
    S2 = S1;
    S3 = S1;
    
    for i=1:nT
        sigma = [sxx(i), sxy(i), sxz(i); sxy(i),syy(i), syz(i); sxz(i), syz(i), szz(i)];
        [rot, sD] = eig(sigma);
        
		sD = diag(sD);
		
		[~, sortI] = sort(sD, 'descend');
		sD = sD(sortI);
		rot = rot(:, sortI);

		if(sum(abs(sD))==0) 
			sD = [1 1 1]';
		end

		relError13 = abs(sD(1)-sD(3))/max(abs(sD));
		relError12 = abs(sD(1)-sD(2))/max(abs(sD));

		if(relError12 < tol)
			sD = sD([3,1,2]);
		end

		if(relError13 < tol) 
			   sD = [1 1 1]';
		else
			sD = 9.*(sD - sD(3))./abs(sD(1) - sD(3)) + 1;
		end
        

        sigma_sqrt = diag(sign(sD).*sqrt(abs(sD)))*rot';

        S1(i, :) = sigma_sqrt(1, :);
        S2(i, :) = sigma_sqrt(2, :);
        S3(i, :) = sigma_sqrt(3, :);       
    end

    function [cost, g] = accuracy(x)
        wCostV = reshape(x, size(w));

        %move w's onto element centers
        wCost = (1/4)*(wCostV(T(:,1),:)+wCostV(T(:,2),:)+wCostV(T(:,3),:)+wCostV(T(:,4),:));

        % rebuild frames and compute matching cost
        % vectorized version 
        expMVec = rodriguesVec(wCost);
		cost = sum( vecNorm(matVecProd([S1 S2 S3],expMVec(:,[2,5,8])))+...
					vecNorm(matVecProd([S1 S2 S3],expMVec(:,[3,6,9]))));
        
        % gradient calculation
        if nargout > 1 % gradient required
            %----- Vectorized version starts here ------ %
            frame2 = matProdVec([S1 S2 S3], expMVec);
            
            %vector components
            vAll2 = bsxfun(@rdivide, frame2(:, [2 5 8]), vecNorm(frame2(:, [2 5 8])));
            vAll3 = bsxfun(@rdivide, frame2(:, [3 6 9]), vecNorm(frame2(:, [3 6 9])));
            
            wCost(vecNorm(wCost) < eps,:) = rand(sum(vecNorm(wCost) < eps),3).*eps; 
            
            thetaAll = vecNorm(wCost);
            thetaiAll = 1./thetaAll;
            
            wmThetaiAll = bsxfun(@times, wCost, thetaiAll);
            ctAll = cos(thetaAll);
            stAll = sin(thetaAll);
            
            KAll = cpMatrixVec(wmThetaiAll);
            KAll2 = matProdSkewSymVec(KAll, KAll);
            
                
            K1All = bsxfun(@times, repmat(thetaiAll,1,6), cpMatrixVec([1 0 0]))...
                    -bsxfun(@times, repmat(wmThetaiAll(:,1),1,6), repmat(thetaiAll,1,6))...
                    .*cpMatrixVec(wmThetaiAll);
            K2All = bsxfun(@times, repmat(thetaiAll,1,6), cpMatrixVec([0 1 0]))...
                    -bsxfun(@times, repmat(wmThetaiAll(:,2),1,6), repmat(thetaiAll,1,6))...
                    .*cpMatrixVec(wmThetaiAll);
            
            K3All = bsxfun(@times, repmat(thetaiAll,1,6), cpMatrixVec([0 0 1]))...
                    -bsxfun(@times, repmat(wmThetaiAll(:,3),1,6), repmat(thetaiAll,1,6))...
                    .*cpMatrixVec(wmThetaiAll);
                
            R1All = skewToFull(repmat(ctAll.*wmThetaiAll(:,1),1, 6).*KAll +...
                  repmat(stAll,1,6).*K1All) + repmat(stAll.*wmThetaiAll(:,1),1,9).*KAll2+repmat((1-ctAll),1,9).*...
                  (matProdSkewSymVec(K1All,KAll)+matProdSkewSymVec(KAll,K1All));
             
            R2All = skewToFull(repmat(ctAll.*wmThetaiAll(:,2),1, 6).*KAll +...
                  repmat(stAll,1,6).*K2All) + repmat(stAll.*wmThetaiAll(:,2),1,9).*KAll2+repmat((1-ctAll),1,9).*...
                  (matProdSkewSymVec(K2All,KAll)+matProdSkewSymVec(KAll,K2All));
              
            R3All = skewToFull(repmat(ctAll.*wmThetaiAll(:,3),1, 6).*KAll +...
                  repmat(stAll,1,6).*K3All) + repmat(stAll.*wmThetaiAll(:,3),1,9).*KAll2+repmat((1-ctAll),1,9).*...
                  (matProdSkewSymVec(K3All,KAll)+matProdSkewSymVec(KAll,K3All));
            
            R2mAllT = [R1All(:,2) R1All(:,5) R1All(:, 8) ...
                      R2All(:,2) R2All(:,5) R2All(:, 8) ...
                      R3All(:,2) R3All(:,5) R3All(:, 8)];
                  
            R3mAllT = [R1All(:,3) R1All(:,6) R1All(:, 9) ...
                      R2All(:,3) R2All(:,6) R2All(:, 9) ...
                      R3All(:,3) R3All(:,6) R3All(:, 9)];
            
            ST = [S1(:,1) S2(:,1) S3(:,1) ...
                  S1(:,2) S2(:,2) S3(:,2) ...
                  S1(:,3) S2(:,3) S3(:,3)];
              
			gVAll = matVecProd(matProdVec(R2mAllT, ST), vAll2) + ...
					matVecProd(matProdVec(R3mAllT, ST), vAll3);
            
            g = 0.25.*V2F'*gVAll;
            g = g(:);
        end
    end

    function [cost, g] = smoothness(x)
        cost = 0.5*(x'*H*x) + 0.5*(x')*x;

        if nargout > 1 % gradient required
            g = H*x + x;
        end
    end


	function [cost,g] = costFunc(x, alpha)
        [B1,B2] = accuracy(x);
        [C1, C2] = smoothness(x);

        if nargout > 1 % gradient required
			g = B2+alpha.*C2;
        end
         
        cost = B1+alpha.*C1;
    end

 
    options = optimoptions(...
        'fmincon',...
        'Display', 'iter',...
        'MaxIterations', 100000,...
        'Algorithm', 'interior-point',...
        'SpecifyObjectiveGradient', true,...
        'Hessian', 'lbfgs',...
        'CheckGradients', false,...
        'OptimalityTolerance', 1e-2);
    
    wOpt = w(:);
    alpha = 10*nT;

    numIter = 30;
    
    trackCost = zeros(numIter,0);

    for iter=1:(numIter)
        disp(iter);
        
        wOpt = fmincon(@(x) costFunc(x, alpha), wOpt(:), [], [], [], [], [], [], [], options);
        trackCost(iter) = accuracy(wOpt);
        alpha = alpha*0.67;
    end
    
    wOpt = reshape(wOpt, size(w));
    wOpt = (1/4)*(wOpt(T(:,1),:)+wOpt(T(:,2),:)+wOpt(T(:,3),:)+wOpt(T(:,4),:));

    dataOut.V = V;
    dataOut.T = T;
    dataOut.stress = stress;
    dataOut.v = zeros(size(wOpt,1),3);
    dataOut.w = zeros(size(wOpt,1),3);
    dataOut.t = zeros(size(wOpt,1),3);

    for ii=1:size(wOpt,1)
        frame = expm([0 -wOpt(ii,3) wOpt(ii,2); wOpt(ii,3) 0 -wOpt(ii,1); -wOpt(ii,2) wOpt(ii,1) 0 ]);
        dataOut.v(ii,:) = frame(1:3,1)';
        dataOut.w(ii,:) = frame(1:3,2)';
        dataOut.t(ii,:) = frame(1:3,3)';
    end

    figure;
    plot(1:numIter,trackCost, 'r');

end

