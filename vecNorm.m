% Vectorized vector 2-norm
function d = vecNorm(V)
	d = sqrt(V(:,1).*V(:,1)+V(:,2).*V(:,2)+V(:,3).*V(:,3));
end

