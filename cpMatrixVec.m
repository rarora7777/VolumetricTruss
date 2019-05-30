% Rotation vector to cross-product matrix (vectorized) conversion
function M = cpMatrixVec(k)
%Output: M: [nx6] = row major matrix storage for n cross-product matrices
%Input:  k: [nx3] = n axes of rotation vectors

%M = [12 13 21 23 31 33]
	M =[-k(:,3) k(:,2) k(:,3) -k(:,1) -k(:,2) k(:,1)];
end