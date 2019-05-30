function plotTriVectors3D(V, F, v)

% Plot stress directly on faces, or interpolate onto vertices and then plot
plotOnVertices = false;
%interpolate onto vertex positions
Xtri = V(:,1);
Ytri = V(:,2);
Ztri = V(:,3);
Xtri2 = Xtri(F);
Ytri2 = Ytri(F);
Ztri2 = Ztri(F);

nF = size(F,1);
CX = sum(Xtri2,2)/3;
CY = sum(Ytri2,2)/3;
CZ = sum(Ztri2,2)/3;

quiver3(CX, CY, CZ, v(:,1), v(:,2), v(:,3));
axis equal
axis vis3d
xlabel('x'); ylabel('y'); zlabel('z')

end

