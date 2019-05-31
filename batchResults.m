%% Batch-process all the results shown in the paper, or select a subset


% Select a subset of problems to solve (see processCSVWithCases).
% Use cases = 1:20 to solve for all the examples in the paper.
% All examples require gptoolbox [Jacobson et al 2018].
% See https://github.com/alecjacobson/gptoolbox
% In addition, cases# 3, 8, 15, and 20 require gptoolbox with tetgen. See
% https://github.com/alecjacobson/gptoolbox/blob/master/wrappers/path_to_tetgen.m
cases = 1:20;

n = numel(cases);

% Specify the value of the resolution parameter `rho` for each problem
resArray = [...
	48, 24, 32, 32, 24,...
	32, 32, 32, 24, 48,...
	32, 32, 32, 24, 32,...
	32, 32, 24, 48, 32];


% Generate all results with beta=0.1 (prefer orthogonality)
resultsOrtho = generateExamples(cases, 0.1*ones(1,n), resArray(cases));

% Generate all results with beta=1.0 (prefer regular spacing)
resultsSpaced = generateExamples(cases, ones(1,n), resArray(cases));
