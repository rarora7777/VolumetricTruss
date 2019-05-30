function [dataOut, timings, mesh] = processCSVwithSim3DCases(example, BETA, res)
    
    useLoubignac = true;
    timings = zeros(1, 4);
    mesh = zeros(1, 2);
    dirichletOnly = false;
       
    % Load geometry and setup loading conditions for the case #`example`.
	% Geometry for cantilever beam and bars is created on the fly.
    % Units are m, kg, N (makes sure meshes/geometries are sane sizes)
    switch example
		case 1 %curved bridge
            [V, T] = readMESH('./data/curved_bridge.mesh');
            
            F = boundary_faces(T);
            deg = accumarray(F(:), 1);
            deg = [deg; zeros(size(V, 1) - numel(deg), 1)];
            
            fixedV = deg > 0 & V(:, 2) - min(V(:, 2)) < 1e-2;
            P = fixedBC(fixedV);
            
            N = normalizerow(cross(V(F(:, 2), :) - V(F(:, 1), :), V(F(:, 3), :) - V(F(:, 2), :), 2));
            f2v = sparse(F(:), repmat((1:size(F, 1))', 3, 1), 1, size(V, 1), size(F, 1));
            NV = (f2v*N)./deg;
            Vb = deg > 0 & NV(:, 2) > max(abs(NV(:, 1)), abs(NV(:, 3))) & NV(:, 2) > 0;
            
            f = singleLoad(size(V, 1), Vb, [0, -1e5./sum(Vb), 0]');
		
		case 2 % mars lander (upper leg)
            [V, T] = readMESH('./data/mars_lander_upper_leg.mesh');
            
            F = boundary_faces(T);
            deg = accumarray(F(:), 1);
            deg = [deg; zeros(size(V, 1) - numel(deg), 1)];
            
            fixedV = deg > 0 & V(:, 1) + V(:, 3) > 9.35;
            P = fixedBC(fixedV);
            
            Vb = deg > 0 & V(:, 1) + V(:, 3) < 5.2;
            f = singleLoad(size(V, 1), Vb, [0, -1e5./sum(Vb), 0]');
            
		case 3 % satellite antenna arm
            [Vobj, Fobj] = readOBJ('./data/antenna.obj');
            [V, T] = tetgen(Vobj, Fobj, 'Flags', sprintf('-q1.2a%0.17f', 0.008*avgedge(Vobj, Fobj)^3/sqrt(2)));
            
            F = boundary_faces(T);
            deg = accumarray(F(:), 1);
            deg = [deg; zeros(size(V, 1) - numel(deg), 1)];
            
            fixedV = deg > 0 & V(:, 3) > 0.04 & max(V(:, 1)) - abs(V(:, 1)) > 1e-4;
            P = fixedBC(fixedV);
            
            Vb = deg > 0 & V(:, 2) > 0.12;
            f = singleLoad(size(V, 1), Vb, [0, -1000./sum(Vb), -500./sum(Vb)]');
		
		case 4 % holey pillar
            [V, T] = readMESH('./data/holey_pillar.mesh');
            
            F = boundary_faces(T);
            deg = accumarray(F(:), 1);
            deg = [deg; zeros(size(V, 1) - numel(deg), 1)];
            
            fixedV = deg > 0 & abs(V(:, 2) - min(V(:, 2))) < 1e-3;
            P = fixedBC(fixedV);
            
            Vb = deg > 0 & (max(V(:, 2)) - V(:, 2)) < 0.1;
            f = singleLoad(size(V, 1), Vb, [0, -500./sum(Vb), 0]');
        
		case 5 % cantilever beam
            %geometry
            dx = linspace(-.1, .1, 24);
            dy = linspace(-0.05, 0.05, 12);
            dz = linspace(-0.05, 0.05, 12);
            [x,y,z] = meshgrid(dx,dy,dz); % a cube
            x = [x(:);0];
            y = [y(:);0];
            z = [z(:);0];
            DT = delaunayTriangulation(x,y,z);
            V = DT.Points;
            T = DT.ConnectivityList; 

            %boundary conditions
            leftV = (V(:,1) == min(V(:,1)));
            rightV = (V(:,1) == max(V(:,1)));
            P = fixedBC(leftV);
            f = singleLoad(size(V,1), rightV, [0, -4 0]');
        
		case 6 % bar under torsion
            %geometry
            dx = linspace(-.1, .1, 24);
            dy = linspace(-0.05, 0.05, 12);
            dz = linspace(-0.05, 0.05, 12);
            [x,y,z] = meshgrid(dx,dy,dz); % a cube
            x = [x(:);0];
            y = [y(:);0];
            z = [z(:);0];
            DT = delaunayTriangulation(x,y,z);
            V = DT.Points;
            T = DT.ConnectivityList; 

            %boundary conditions
            leftV = (V(:,1) == min(V(:,1)));
            rightV = (V(:,1) == max(V(:,1)));
            Vb = V(:, 2) == max(V(:, 2));
            fixedV = V(:, 2) == min(V(:, 2));
            frontV = V(:, 3) == max(V(:, 3));
            backV = V(:, 3) == max(V(:, 3));
            
            P = fixedBC(leftV);
            
            mult = 4.0;
            f1 = singleLoad(size(V,1), ...
                (rightV & backV),...
                [0 mult 0]');
            f2 = singleLoad(size(V,1), ...
                (rightV & frontV),...
                [0 -mult 0]');
            f3 = singleLoad(size(V,1), ...
                (rightV & Vb) & ~(backV | frontV),...
                [0 0 mult]');
            f4 = singleLoad(size(V,1), ...
                (rightV & fixedV) & ~(backV | frontV),...
                [0 0 -mult]');
            
            f = f1 + f2 + f3 + f4;
        
		case 7 % bar under tension (pulled from both ends)
            dx = linspace(-.1, .1, 24);
            dy = linspace(-0.05, 0.05, 12);
            dz = linspace(-0.05, 0.05, 12);
            [x,y,z] = meshgrid(dx,dy,dz); % a cube
            x = [x(:);0];
            y = [y(:);0];
            z = [z(:);0];
            DT = delaunayTriangulation(x,y,z);
            V = DT.Points;
            T = DT.ConnectivityList; 

            %boundary conditions
            leftV = (V(:,1) == min(V(:,1)));
            rightV = (V(:,1) == max(V(:,1)));
            
            P = fixedBC(leftV);
            
            mult = 40.0;
            f = singleLoad(size(V,1), ...
                rightV,...
                [mult 0 0]');
			
		case 8 % simple bridge
            %load data
            [Vobj,Fobj] = readOBJ('./data/bridge.obj');
            [V,T] = tetgen(Vobj,Fobj,'Flags',sprintf('-q1.2a%0.17f',16*avgedge(Vobj,Fobj)^3/(6*sqrt(2))));
            
            %fix bottom of bridge
            Vb = (V(:,2) == max(V(:,2)));
            fixedV = (V(:,2) == min(V(:,2)));
            P = fixedBC(fixedV);
            
            f = singleLoad(size(V,1), Vb, [0, -1e4./sum(Vb) 0]');
        
		case 9 % mars lander (lower leg)
            [V, T] = readMESH('./data/mars_lander_lower_leg.mesh');
            
            F = boundary_faces(T);
            deg = accumarray(F(:), 1);
            deg = [deg; zeros(size(V, 1) - numel(deg), 1)];
            
            fixedV = deg > 0 & V(:, 2) < -2.8;
            P = fixedBC(fixedV);
            
            Vb = deg > 0 & V(:, 2) > -0.49;
            f = singleLoad(size(V, 1), Vb, [0, -1e5./sum(Vb), 0]');
		
		case 10 % pavilion
            [V, T] = readMESH('./data/pavilion.mesh');
            
            F = boundary_faces(T);
            deg = accumarray(F(:), 1);
            deg = [deg; zeros(size(V, 1) - numel(deg), 1)];
            
            fixedV = deg > 0 & abs(V(:, 2) - min(V(:, 2))) < 1e-3;
            P = fixedBC(fixedV);
            
            Vb = deg > 0 & (max(V(:, 2)) - V(:, 2)) < 1e-3;
            f = singleLoad(size(V, 1), Vb, [0, -1e5./sum(Vb), 0]');
        
		case 11 % bookcase
            [V, T] = readMESH('./data/bookcase.mesh');
            
            F = boundary_faces(T);
            deg = accumarray(F(:), 1);
            deg = [deg; zeros(size(V, 1) - numel(deg), 1)];
            
            N = normalizerow(cross(V(F(:, 2), :) - V(F(:, 1), :), V(F(:, 3), :) - V(F(:, 2), :), 2));
            f2v = sparse(F(:), repmat((1:size(F, 1))', 3, 1), 1, size(V, 1), size(F, 1));
            NV = (f2v*N)./deg;
            
            fixedV = deg > 0 & abs(V(:, 3) - min(V(:, 3))) < 1e-3;
            P = fixedBC(fixedV);
            
            Vb1 = deg > 0 & V(:, 3) > 0 & V(:, 2) < 0.4 & NV(:, 2) > 0 & NV(:, 2) > max(abs(NV(:, 1)), abs(NV(:, 3)));
            f1 = singleLoad(size(V, 1), Vb1, [0, -5e2./sum(Vb1), 0]');
            
            Vb2 = deg > 0 & V(:, 3) > 0 & V(:, 1) < 0.01 & V(:, 1) > -0.1 & V(:, 2) < 0.45 & V(:, 2) > 0.2;
            f2 = singleLoad(size(V, 1), Vb2, [1e2./sum(Vb2), 0, 0]');
            
            Vb3 = deg > 0 & V(:, 3) > 0 & V(:, 1) < -0.3 & V(:, 1) > -0.37 & V(:, 2) < 0.11 & V(:, 2) > -0.12;
            f3 = singleLoad(size(V, 1), Vb3, [-1e2./sum(Vb3), 0, 0]');
            
            Vb4 = deg > 0 & V(:, 3) > 0 & V(:, 1) < 0.3 & V(:, 1) > 0.2 & V(:, 2) < -0.15 & V(:, 2) > -0.42;
            f4 = singleLoad(size(V, 1), Vb4, [-1e2./sum(Vb4), 0, 0]');
            
            f = f1 + f2 + f3 + f4;
            
        case 12 % mars lander (body)
            [V, T] = readMESH('./data/mars_lander_body.mesh');
            
            F = boundary_faces(T);
            deg = accumarray(F(:), 1);
            deg = [deg; zeros(size(V, 1) - numel(deg), 1)];
            
            midX = (max(V(:, 1)) + min(V(:, 1)))/2;
            midZ = (max(V(:, 3)) + min(V(:, 3)))/2;
            
            fixedV = deg > 0 & abs(V(:, 1) - midX) + abs(V(:, 3) - midZ) > 3.73;
            P = fixedBC(fixedV);
            
            Vb = deg > 0 & abs(V(:, 1)) < 0.82 & abs(V(:, 2) + 0.8) < 0.02;
            f = singleLoad(size(V, 1), Vb, [0, -1e5./sum(Vb), 0]');

		case 13 % arched bridge
            [V, T] = readMESH('./data/arched_bridge.mesh');
            
            F = boundary_faces(T);
            deg = accumarray(F(:), 1);
            deg = [deg; zeros(size(V, 1) - numel(deg), 1)];
            
            fixedV = deg > 0 & V(:, 2) - min(V(:, 2)) < 0.05;
            P = fixedBC(fixedV);
            
            Vb = deg > 0 & (3*V(:, 1) - 5*V(:, 2) < -17 | V(:, 1) + 2*V(:, 2) > 6 | V(:, 2) > 2);
            f = singleLoad(size(V, 1), Vb, [0, -1e4./sum(Vb), 0]');
        
		case 14 %quadcopter frame
            %load data
            [V, T] = readMESH('./data/quadcopter_frame.mesh');
            
            % fix top edge (assume propellors are like a fixed boundary
            % condition)
            fixedV = (abs(V(:,2)- max(V(:,2))) < 0.001);
            Vb = (abs(V(:,2)- min(V(:,2))) < 0.001);
            P = fixedBC(fixedV);
            f = singleLoad(size(V,1), Vb, [0, -100./sum(Vb) 0]');
		
		case 15 % helicopter top pylon
            [Vobj, Fobj] = readOBJ('./data/top_pylon.obj');
            [V, T] = tetgen(Vobj, Fobj, 'Flags', sprintf('-q1.2a%0.17f', 8*avgedge(Vobj, Fobj)^3/sqrt(2)));
            
            F = boundary_faces(T);
            deg = accumarray(F(:), 1);
            deg = [deg; zeros(size(V, 1) - numel(deg), 1)];
            
            fixedV = deg > 0 & V(:, 1) > 2.38 & V(:, 1) < 3.1 & V(:, 2) > .1 & abs(V(:, 3)) < .25;
            P = fixedBC(fixedV);
            
            Vb = deg > 0 & V(:, 2) < .03;
            f = singleLoad(size(V, 1), Vb, [0, -1e4./sum(Vb), 0]');
        
		case 16 % chair under sitting load
            [V, T] = readMESH('./data/chair.mesh');
            
            F = boundary_faces(T);
            deg = accumarray(F(:), 1);
            deg = [deg; zeros(size(V, 1) - numel(deg), 1)];
            
            fixedV = deg > 0 & abs(V(:, 2) - min(V(:, 2))) < 1e-3*(max(V(:, 2)) - min(V(:, 2)));
            P = fixedBC(fixedV);
            
            Vb = deg > 0 & V(:, 2) > 0.05 & V(:, 2) < 0.06 & V(:, 3) > -0.12;
            f = singleLoad(size(V, 1), Vb, [0, -700./sum(Vb), 0]');
        
		case 17 % chair with rocking load
            [V, T] = readMESH('./data/chair.mesh');
            
            F = boundary_faces(T);
            deg = accumarray(F(:), 1);
            deg = [deg; zeros(size(V, 1) - numel(deg), 1)];
            
            fixedV = deg > 0 & abs(V(:, 2) - min(V(:, 2))) < 1e-3*(max(V(:, 2)) - min(V(:, 2)));
            P = fixedBC(fixedV);
            
            Vb1 = deg > 0 & V(:, 2) > 0.05 & V(:, 2) < 0.06 & V(:, 3) > -0.12;
            f1 = singleLoad(size(V, 1), Vb1, [0, -500./sum(Vb1), 0]');
            
            Vb2 = deg > 0 & V(:, 2) > 0.06 & V(:, 3) < -0.11 & V(:, 3) > -0.12;
            f2 = singleLoad(size(V, 1), Vb2, [0, 0, -500./sum(Vb2)]');
            
            f = f1 + f2;
			
		case 18 %climbing_hold
            [V, T] = readMESH('./data/climbing_hold.mesh');
            
            F = boundary_faces(T);
            deg = accumarray(F(:), 1);
            deg = [deg; zeros(size(V, 1) - numel(deg), 1)];
            
            fixedV = deg > 0 & V(:, 3) - min(V(:, 3)) < 1e-3;
            P = fixedBC(fixedV);
            
            Vb = deg > 0 & V(:, 2) < -.10 & V(:, 2) > -.40 & 2*V(:, 3) - V(:, 2) < .53;
            
            f = singleLoad(size(V, 1), Vb, [0, -700./sum(Vb), 0]');
		
		case 19 % holey sculpture
            [V, T] = readMESH('./data/holey_sculpture.mesh');
            
            F = boundary_faces(T);
            deg = accumarray(F(:), 1);
            deg = [deg; zeros(size(V, 1) - numel(deg), 1)];
            
            fixedV = deg > 0 & V(:, 2) - min(V(:, 2)) < 0.01;
            
            P = fixedBC(fixedV);
            
            Vb = deg > 0 & ~fixedV;
            
            f = singleLoad(size(V, 1), Vb, [1e5./sum(Vb), 0, 0]');
            
        case 20 % jet engine bracket
            [Vobj, Fobj] = readOBJ('./data/engine_bracket.obj');
            [V, T] = tetgen(Vobj, Fobj, 'Flags', sprintf('-q1.2a%0.17f', 8*avgedge(Vobj, Fobj)^3/sqrt(2)));
            
            F = boundary_faces(T);
            deg = accumarray(F(:), 1);
            deg = [deg; zeros(size(V, 1) - numel(deg), 1)];
            
            fixedV = deg > 0 & V(:, 2) > -3.5 & V(:, 2) < -2 &...
                ((V(:, 1) > -10.6 & V(:, 1) < -8 & V(:, 3) > -3 & V(:, 3) < 6) |...
                (V(:, 1) > 8 & V(:, 1) < 10.6 & V(:, 3) > -3.5 & V(:, 3) < 4));
            P = fixedBC(fixedV);
            
            Vb = deg > 0 & (V(:, 2)-1.7).^2 + (V(:, 3)+4.55).^2 < 1.25*1.25;
            f = singleLoad(size(V, 1), Vb, [0, (4.3*1e4/sqrt(2))./sum(Vb), -(4.3*1e4/sqrt(2))./sum(Vb)]');
        
		% Add your own problems here
		% case 21 % custom problem
		
        otherwise
            dataOut = [];
            timings = [];
            mesh = [];
            return
    end
    
	
%% Compute the stress field using GAUSS [Levin et al 2017]
    mesh(1) = size(V,1);
    mesh(2) = size(T,1);
    
	tic;
    fem = WorldFEM('elastic_linear_tetrahedra', V, T);
    K = stiffness(fem);
    
    %solve static problem
    Kbc = P*K*P';
    
    if ~dirichletOnly
        f = f + force(fem);
        u = P'*(Kbc\(P*f));
    end
    
    

	% Loubignac iterations for smoothing the stress field
    if useLoubignac
        Stmp = loubignac(fem, P, f, u, 1e-4, 1000);
        data.vertexStress = Stmp;

        %Average Stress back onto elements
        ti = (1:size(T,1))';
        bary = repmat(0.25, size(T,1),4);

        Stmp1 = Stmp(:,1);
        Stmp1 = Stmp1(T(ti,:));
        Stmp1 = dot(bary',Stmp1')';

        Stmp2 = Stmp(:,2);
        Stmp2 = Stmp2(T(ti,:));
        Stmp2 = dot(bary',Stmp2')';

        Stmp3 = Stmp(:,3);
        Stmp3 = Stmp3(T(ti,:));
        Stmp3 = dot(bary',Stmp3')';

        Stmp4 = Stmp(:,4);
        Stmp4 = Stmp4(T(ti,:));
        Stmp4 = dot(bary',Stmp4')';

        Stmp5 = Stmp(:,5);
        Stmp5 = Stmp5(T(ti,:));
        Stmp5 = dot(bary',Stmp5')';

        Stmp6 = Stmp(:,6);
        Stmp6 = Stmp6(T(ti,:));
        Stmp6 = dot(bary',Stmp6')';

        Stmp = [Stmp1 Stmp2 Stmp3 Stmp4 Stmp5 Stmp6];
    else
	% Standard linear FEM (without Loubignac)
        Stmp = stress(fem, u);
    end
  
    timings(1) = toc;

    data.V = V;
    data.T = T;
    % our stress format is (xx, yy, zz, xy, yz, xz)
    data.stress = [Stmp(:,1); Stmp(:,2); Stmp(:,3); Stmp(:,4); Stmp(:,5); Stmp(:,6)];
    
%% Plot the stress field (optional)
    figure; plotStressEigs3D(data);
    drawnow;
	
%% Solve for smooth stress-aligned frame field
    tic;
    dataFrames = fitFramesToData3D(data.V, data.T, data.stress);
    timings(2) = toc;

%% Solve for texture parametrization
    tic;
    dataTex = fitTexCoords3D(dataFrames, BETA);
    timings(3) = toc;
    
    dataOut = dataTex;
    dataOut.u = matrixnormalize(dataOut.u);

%% Extract the truss layout
	
	% default resolution is 32
    if nargin < 3
        res = 32;
    end
	
    tic;
	dataOut = tex2CurvesTet(dataOut, res, false);
    timings(4) = toc;
   
    clear fem;
end