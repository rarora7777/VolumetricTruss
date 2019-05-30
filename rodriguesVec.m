function M = rodriguesVec(W)
%vecrorized implementation of the rodrigues rotation formula
%input W - nx3 martrix of angular rotation vectors [wx, wy, wz]
%output M -nx9 list of ourput matrix coefficients in row major order

%compute magnitudes
mag = sqrt(W(:,1).*W(:,1)+W(:,2).*W(:,2)+W(:,3).*W(:,3));
%compute normalized W's
W = bsxfun(@rdivide, W, mag);
W(isnan(W)) = 0;
st = sin(mag);
ct = cos(mag);
K = cpMatrixVec(W);
K2 = bsxfun(@times, matProdSkewSymVec(K,K), (1-ct));
K = bsxfun(@times, K,st); 
M = [1+K2(:,1), K(:,1)+K2(:,2), K(:,2)+K2(:,3), K(:,3)+K2(:,4), 1+K2(:,5), K(:,4)+K2(:,6), K(:,5)+K2(:,7), K(:,6)+K2(:,8), 1+K2(:,9)];

end

