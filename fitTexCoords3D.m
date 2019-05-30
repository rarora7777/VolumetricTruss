% this function takes, as input, the output of the vector field fitting phase
% and parameterizes it by solving for a smooth, coordinate aligned piecewise
% linear map 
function dataOut = fitTexCoords3D(dataFit, BETA)
	V = dataFit.V;
	T = dataFit.T;

	%1. Initialize Variables
	v = dataFit.v;
	w = dataFit.w;
	t = dataFit.t;

	%Gradient Operators
	nV = size(V, 1);
	nT = size(T, 1);

	G = grad(V, T);
	Gx = G(1:nT, 1:nV);
	Gy = G(nT+1:2*nT, 1:nV);
	Gz = G(2*nT+1:3*nT, 1:nV);

	vG = sparse(bsxfun(@times, v(:,1), Gx) + bsxfun(@times, v(:,2), Gy) + bsxfun(@times, v(:,3), Gz));
	wG = sparse(bsxfun(@times, w(:,1), Gx) + bsxfun(@times, w(:,2), Gy) + bsxfun(@times, w(:,3), Gz));
	tG = sparse(bsxfun(@times, t(:,1), Gx) + bsxfun(@times, t(:,2), Gy) + bsxfun(@times, t(:,3), Gz));


	zb = sparse(size(vG,1), size(vG,2));	% block of zeros
	zv = zeros(size(vG,1),1);	% vector of zeros
	ov = ones(size(vG,1),1);	% vector of ones


	% First term of Eq. (11): Equal-spacing
	Aopt = sparse([...
		vG zb zb;...
		zb wG zb;...
		zb zb tG...
	]);
	
	bopt = [ov; ov; ov];

	% Second term of Eq. (11): Alignment
	Aeq = sparse([...
		zb vG zb;...
		zb zb vG;...
		wG zb zb;...
		zb zb wG;...
		tG zb zb;...
		zb tG zb;...
		1 zb(1, 1:end-1) zb(1,:) zb(1,:); zb(1,:) 1 zb(1,1:end-1) zb(1,:);...
		zb(1,:) zb(1,:) 1 zb(1,1:end-1)...
	]);
	
	beq = [zv; zv; zv; zv; zv; zv; 0; 0; 0];
	
	% Solve for the parametrization
	param = quadprog(BETA*(Aopt'*Aopt) + Aeq'*Aeq, -BETA.*Aopt'*bopt - Aeq'*beq, [],[], [],[]);

	% Save results and return
	dataOut = dataFit;
	dataOut.V = V;
	dataOut.T = T;
	dataOut.u = [param(1:nV) param((1+nV):2*nV) param((1+2*nV):3*nV)];
	dataOut.u = matrixnormalize(dataOut.u);
	dataOut.Beta = BETA;
end
