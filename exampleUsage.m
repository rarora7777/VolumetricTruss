exampleIdx = 5; % cantilever beam
resolution = 16; % select `rho=16`
b = 0.1;  % select `beta=0.1` for parametrization solve

truss = processCSVwithSim3DCases(exampleIdx, b, resolution);

% Generate sparser structure
n = 2;  % generate sparse structure by selecting every n-th curve for each parameter
sparseTruss = sparsifyElements(truss, n);

% Visualize the sparse truss layout
drawTrussGraph(sparseTruss);

% Visualize cylinder-extruded truss and save generated meshes to OBJ files
rE = 2e-3;  % radius of cylindrical elements: 2mm
rN = 3e-3;  % radius of spherical nodes: 3mm
generateNodeElemMeshes(sparseTruss, './results/beam', rE, rN);

% Generate manifold mesh using libigl's booleans

% first, collapse tiny elements since the wire_mesh code does not tolerate
% near-degenerate elements
d = collapseSmallEdges(sparseTruss, 1e-3);
node = d.Node;
elem = d.Elem;
% also, remove any disconnected nodes, possibly formed during the last
% operation
[node, SVI] = remove_unreferenced(node, elem);
elem = SVI(elem);
% Generate solid mesh using libigl/gptoolbox booleans
[VT, FT] = wire_mesh(node, elem, 'PolySize', 6, 'Thickness', 3e-3, 'Solid', true);
% Visualize and save the mesh
tsurf(FT, VT, 'EdgeColor', 'k', 'FaceColor', [203, 109, 81]/255, 'EdgeAlpha', 0.5);
axis equal off;

writeOBJ('./results/beam_solid.obj', VT, FT);

clear node elem d SVI;