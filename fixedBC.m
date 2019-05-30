function P = fixedBC(bc)
    n = 1:numel(bc);
    P = sparse(kron(sparse(n,n, abs(1-bc)), sparse(eye(3,3))));
    P(sum(P,2) == 0, :) = [];
end