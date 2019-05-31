function M = skewToFull(S)
zz = zeros(size(S,1),1);
M = [zz S(:,1) S(:,2) S(:,3) zz S(:,4) S(:,5) S(:,6) zz];
end

