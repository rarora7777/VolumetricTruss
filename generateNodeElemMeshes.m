function generateNodeElemMeshes(data, outPath, rE, rN)
% Generate and save truss meshes suitable for rendering
% Input
% data: struct with fields Node and Elem
% outPath: path (including filename) for saving the output meshes
% rE: radius of the cylindrical elements
% rN: radius of the spherical nodes
%
% Output
% Visualizes the mesh and saves two OBJ files:
% <outPath>_nodes.obj, <outPath>_elem.obj
%
% WARNING: These meshes are only suitable for rendering.
% Use gptoolbox's wire_mesh function to generate a boolean'd manifold mesh

    nodePath = [outPath, '_nodes.obj'];
    elemPath = [outPath, '_elems.obj'];
    
    node = data.Node;
    elem = data.Elem;
    
    
    [sV, sF] = readOBJ('./util_meshes/sphere_8_5.obj');
    [cV, cF] = readOBJ('./util_meshes/cylinder_6.obj');

    nsV = size(sV, 1);
    nsF = size(sF, 1);
    ncV = size(cV, 1);
    ncF = size(cF, 1);

    sV = [sV, ones(nsV, 1)];
    cV = [cV, ones(ncV, 1)];
    
    [node, SVI, ~] = remove_unreferenced(node, elem);
    elem = SVI(elem);


    n = size(node, 1);
    m = size(elem, 1);

    nodeF = zeros(nsF*n, 3);
    nodeV = zeros(nsV*n, 4);

    elemF = zeros(ncF*m, 3);
    elemV = zeros(ncV*m, 4);

    for i=1:n
		S = diag([rN rN rN 1])';
        T = [[eye(3); ones(1, 3)], [node(i, :)'; 1]]';
        nodeV((i-1)*nsV + (1:nsV), :) = sV*S*T;
        nodeF((i-1)*nsF + (1:nsF), :) = sF + (i-1)*nsV;
    end

    eVec = node(elem(:, 1), :) - node(elem(:, 2), :);

    eLen = sqrt(sum(eVec.^2, 2));
    eVec = eVec./eLen;

    vec2 = cross(eVec, repmat([0 1 0], m, 1));
    yIdx = sum(vec2.^2, 2) < 1e-6;
    vec2(yIdx, :) = cross(eVec(yIdx, :), repmat([0 0 1], sum(yIdx), 1));
    vec2 = vec2./sqrt(sum(vec2.^2, 2));

    vec3 = cross(eVec, vec2);
    vec3 = vec3./sqrt(sum(vec3.^2, 2));

    ePos = (node(elem(:, 1), :) + node(elem(:, 2), :))/2;

    for i=1:m
		S = diag([eLen(i) rE rE 1])';
        T = [[eye(3); ones(1, 3)], [ePos(i, :)'; 1]]';

        R = [[eVec(i, :), 0]; [vec2(i, :), 0]; [vec3(i, :), 0]; [0 0 0 1]];
        elemV((i-1)*ncV + (1:ncV), :) = cV*S*R*T;
        elemF((i-1)*ncF + (1:ncF), :) = cF + (i-1)*ncV;
    end

    figure
    tsurf(nodeF, nodeV(:, 1:3), 'FaceColor', [203, 109, 81]/255, 'EdgeAlpha', 0);
    hold on
    tsurf(elemF, elemV(:, 1:3), 'FaceColor', [201, 192, 187]/255, 'EdgeAlpha', 0);
    hold off
    axis equal
    light
    
	writeOBJ(nodePath, nodeV(:, 1:3), nodeF);
	writeOBJ(elemPath, elemV(:, 1:3), elemF);
end