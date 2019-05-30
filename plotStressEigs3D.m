%plot stress vector field
function plotStressEigs3D(data, useLinearLength)

% Plot stress directly on faces, or interpolate onto vertices and then plot
plotOnVertices = false;
%1. eigendecomposition of each tensor
v1 = zeros(size(data.V,1), 2);
v2 = zeros(size(data.V, 1), 2);
v3 = zeros(size(data.V, 1), 2);

%interpolate onto vertex positions
Xtri = data.V(:,1);
Ytri = data.V(:,2);
Ztri = data.V(:,3);
Xtri2 = Xtri(data.T);
Ytri2 = Ytri(data.T);
Ztri2 = Ztri(data.T);

nF = size(data.T,1);
CX = sum(Xtri2,2)/3;
CY = sum(Ytri2,2)/3;
CZ = sum(Ztri2,2)/3;

Xtri = CX;
Ytri = CY;
Ztri = CZ;
sxx = data.stress(1:nF);
syy = data.stress(nF + (1:nF));
szz = data.stress(2*nF + (1:nF));
sxy = data.stress(3*nF + (1:nF));
syz = data.stress(4*nF + (1:nF));
sxz = data.stress(5*nF + (1:nF));

v1 = zeros(size(CX,1), 3);
v2 = zeros(size(CX, 1), 3);
v3 = zeros(size(CZ, 1), 3);

ratio = zeros(size(CX, 1), 2);

for ii=1:numel(Xtri)
    S = [sxx(ii) sxy(ii) sxz(ii); sxy(ii) syy(ii) syz(ii); sxz(ii) syz(ii) szz(ii)];
    S(isnan(S)) = 0;
    [U,D] = eig(S);
    if(sum(abs(D)) < 1e-8) 
        v1(ii,:) = [0 0 0];
        v2(ii, :) = [0 0 0];
        v3(ii, :) = [0 0 0];
    else
        Dvec = diag(D);
        [Y,I] = sort(abs(Dvec), 'descend');
        U = U(:, I);
        if nargin > 1 && useLinearLength
            v1(ii,:) = (Y(1))*U(:, 1)';
            v2(ii, :) = (Y(2))*U(:, 2)';
            v3(ii, :) = (Y(3))*U(:, 3)';
        else
            v1(ii,:) = log(Y(1))*U(:, 1)';
            v2(ii, :) = log(Y(2))*U(:, 2)';
            v3(ii, :) = log(Y(3))*U(:, 3)';
        end
        
        ratio(ii, 1) = real(Y(1))/real(Y(2));
        ratio(ii, 2) = real(Y(2))/real(Y(3));
    end
end

    hold on

    %plot the drawn curve: disabled for now
%     C = data.Curves;
%     for ii = 1:numel(C)
%         line([C{ii}(1:(end-1),1)';C{ii}(2:end,1)'],[C{ii}(1:(end-1),2)';C{ii}(2:end,2)'] )
%     end

quiver3(Xtri, Ytri, Ztri, v1(:,1), v1(:,2), v1(:,3));
quiver3(Xtri, Ytri, Ztri, v2(:,1), v2(:,2), v2(:,3));
quiver3(Xtri, Ytri, Ztri, v3(:,1), v3(:,2), v3(:,3));
axis equal
axis vis3d
hold off

    %store everything in a vtkwrite compatible format
%     SXX =  data.stress(1:nF);
%     SYY = data.stress((nF+1):2*nF);
%     SXY = data.stress((2*nF+1):3*nF);
%     ZERO = zeros(size(SXX,1),1);
%     ONE = ones(size(SXX,1),1);
%     tx = [SXX'; SXY'; ZERO'];
%     ty = [SXY'; SYY'; ZERO'];
%     tz = [ZERO'; ZERO'; ONE'];
% 
%     tx = tx(:);
%     ty = ty(:);
%     tz = tz(:);
end

