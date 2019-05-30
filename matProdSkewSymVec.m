function M = matProdSkewSymVec(A,B) 
%Input: A [nx6] - skew symmetric matrix (missing zeros on the diagonal)
%       B [nx6] - skew symmetric matrix 
	M = [...
		A(:,1).*B(:,3)+A(:,2).*B(:,5),... %11
		A(:,2).*B(:,6),... %12
		A(:,1).*B(:,4),... %13
		A(:,4).*B(:,5),... %21
		A(:,3).*B(:,1)+A(:,4).*B(:,6),... %22
		A(:,3).*B(:,2),... %23
		A(:,6).*B(:,3),... %31
		A(:,5).*B(:,1),... %32
		A(:,5).*B(:,2)+A(:,6).*B(:,4)]; %33
end
