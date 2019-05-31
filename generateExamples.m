function [data, timings, mesh, errors] = generateExamples(cases, beta, res)
	% takes in a list of case ids and resolutions, returns a cell array with all the data
    cases = sort(cases);
    ii = 1;
    data = cell(numel(cases), 1);
    mesh = cell(numel(cases), 1);
    timings = cell(numel(cases), 1);
    errors = cell(numel(cases), 1);
	
    for jj = cases
        try
            [data{ii}, timings{ii}, mesh{ii}] = processCSVwithSim3DCases(jj, beta(ii), res(ii));
            disp([mesh{ii} timings{ii}]);
            result = data{ii};
            time = timings{ii};
            problemSize = mesh{ii};
            save(['./results/res_', num2str(beta(ii)*10), '_', num2str(jj), '_', num2str(res(ii)), '.mat'],...
                'result', 'time', 'problemSize');
            figure; drawTrussGraph(data{ii});
        catch ME
            errors{ii} = ME;
            disp(ME.identifier);
            disp(ME.message);
        end
        ii = ii + 1;
    end
end