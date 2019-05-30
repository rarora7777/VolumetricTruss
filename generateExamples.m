function [data, timings, mesh, errors] = generateExamples(cases, beta, res)
	% takes in a list of case ids and resolutions, returns a cell array with all the data
    cases = sort(cases);
    ii = 0;
    data = cell(numel(cases), 1);
    mesh = cell(numel(cases), 1);
    timings = cell(numel(cases), 1);
    errors = cell(numel(cases), 1);
	
    for jj = cases
        try
            [data{ii+1}, timings{ii+1}, mesh{ii+1}] = processCSVwithSim3DCases(jj, beta(ii+1), res(ii+1));
            disp([mesh{ii+1} timings{ii+1}]);
            result = data{ii+1};
            time = timings{ii+1};
            problemSize = mesh{ii+1};
            save(['./results/res_', num2str(beta(ii+1)*10), '_', num2str(jj), '_', num2str(res(ii+1)), '.mat'],...
                'result', 'time', 'problemSize');
            figure; drawTrussGraph(data{ii+1});
        catch ME
            errors{ii+1} = ME;
            disp(ME.identifier);
            disp(ME.message);
        end
        ii = ii + 1;
    end
end