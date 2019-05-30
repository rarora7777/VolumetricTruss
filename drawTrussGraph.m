function drawTrussGraph(dataIn, width, colour, extraParam)
% Super-fast function for visualizing a truss by treating each element
% as a degenerate triangle
    if isstruct(dataIn)
        node = dataIn.Node;
        elem = dataIn.Elem;
        if nargin == 1
            width = 1;
        end
        if nargin < 3
            colour = 'k';
        end
    else
        node = dataIn;
        elem = width;
        if nargin < 3
            width = 1;
        else
            width = colour;
        end
        if nargin < 4
            colour = 'k';
        else
            colour = extraParam;
        end
    end
    
    if size(node, 2) == 3   % 3D
        tsurf(elem(:, [1 1 2]), node, 'LineWidth', width, 'EdgeColor', colour);
    else    % 2D
        tsurf(elem(:, [1 1 2]), [node, zeros(size(node, 1), 1)], 'LineWidth', width, 'EdgeColor', colour);
    end
    axis equal off
end