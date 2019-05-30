% Vectorized matrix-vector product
function V = matVecProd(A,X)
% Input: A [nx9] - each row is a flattened 3x3 matrix
%        X [nx3] - each row is a 3-vector

	V = [...
		A(:,1).*X(:,1)+A(:,2).*X(:,2)+A(:,3).*X(:,3),...
		A(:,4).*X(:,1)+A(:,5).*X(:,2)+A(:,6).*X(:,3),...
		A(:,7).*X(:,1)+A(:,8).*X(:,2)+A(:,9).*X(:,3)];
end
    

