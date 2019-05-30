function f = singleLoad(nV, bc, pointForce)
    f = zeros(nV,1); 
    f(bc) = 1;
    f = kron(f, pointForce);
end