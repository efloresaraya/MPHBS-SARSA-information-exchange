function out = FDB_TLABC_cec2022_mod(fid, cfg)
% FDB-TLABC adapted to CEC2022 with exact FE budget control.
%
% Usage:
%   out = FDB_TLABC_cec2022_mod(fid, cfg)
%
% Required:
%   fid          : CEC2022 function id
%   cfg.pop      : population size
%   cfg.dim      : problem dimension
%   cfg.maxFE    : exact FE budget
%   cfg.lb       : scalar or 1-by-D lower bound
%   cfg.ub       : scalar or 1-by-D upper bound
%
% Optional:
%   cfg.limit    : scout limit, default 200
%   cfg.CR       : crossover rate, default 0.5
%
% Output:
%   out.best_f   : best raw objective value found
%   out.best_x   : best solution
%   out.curve    : best-so-far raw objective curve
%   out.fe_used  : exact number of function evaluations used

    if nargin < 2
        error('FDB_TLABC_cec2022_mod requires fid and cfg.');
    end
    if isempty(which('cec22_test_func'))
        error('Required function not found on MATLAB path: cec22_test_func');
    end

    cfg = fill_defaults(cfg);

    popsize = cfg.pop;
    D       = cfg.dim;
    maxFES  = cfg.maxFE;

    Xmin = expand_bound(cfg.lb, D);
    Xmax = expand_bound(cfg.ub, D);

    if numel(Xmin) ~= D || numel(Xmax) ~= D
        error('Bounds must be scalar or length-D vectors.');
    end

    obj = @(x) cec22_eval_row(x, fid);

    trial = zeros(1, popsize);
    limit = cfg.limit;
    CR    = cfg.CR;

    X = repmat(Xmin, popsize, 1) + rand(popsize, D) .* repmat(Xmax - Xmin, popsize, 1);
    val_X = zeros(1, popsize);

    FES = 0;
    for i = 1:popsize
        if FES >= maxFES
            break;
        end
        val_X(i) = obj(X(i,:));
        FES = FES + 1;
    end

    [val_gBest, min_index] = min(val_X);
    gBest = X(min_index(1), :);
    curve = val_gBest;

    while FES < maxFES

        % Teaching-based employed bee phase
        for i = 1:popsize
            if FES >= maxFES
                break;
            end

            [~, sortIndex] = sort(val_X);
            mean_result = mean(X, 1);
            Best = X(sortIndex(1), :);
            TF = round(1 + rand*(1));

            Xi = X(i,:) + (Best - TF*mean_result).*rand(1, D);

            % Diversity learning
            r = generateR_local(popsize, i);
            F = rand;
            V = X(r(1),:) + F*(X(r(2),:) - X(r(3),:));
            flag = (rand(1, D) <= CR);
            Xi(flag) = V(flag);

            Xi = boundary_repair_local(Xi, Xmin, Xmax, 'reflect');

            val_Xi = obj(Xi);
            FES = FES + 1;

            if val_Xi < val_X(i)
                val_X(i) = val_Xi;
                X(i,:)   = Xi;
                trial(i) = 0;
            else
                trial(i) = trial(i) + 1;
            end
        end

        % Learning-based onlooker bee phase
        for k = 1:popsize
            if FES >= maxFES
                break;
            end

            Fitness = calculateFitness_local(val_X);
            i = fitnessDistanceBalance_local(X, Fitness);
            j = randi(popsize);
            while j == i
                j = randi(popsize);
            end

            if val_X(i) < val_X(j)
                Xi = X(i,:) + rand(1, D).*(X(i,:) - X(j,:));
            else
                Xi = X(i,:) + rand(1, D).*(X(j,:) - X(i,:));
            end

            Xi = boundary_repair_local(Xi, Xmin, Xmax, 'reflect');

            val_Xi = obj(Xi);
            FES = FES + 1;

            if val_Xi < val_X(i)
                val_X(i) = val_Xi;
                X(i,:)   = Xi;
            end
        end

        % Generalized oppositional scout bee phase
        ind = find(trial == max(trial), 1, 'first');
        if ~isempty(ind) && trial(ind) > limit && FES < maxFES
            trial(ind) = 0;

            sol     = (Xmax - Xmin).*rand(1, D) + Xmin;
            solGOBL = (max(X,[],1) + min(X,[],1)).*rand(1, D) - X(ind,:);
            newSol  = [sol; solGOBL];
            newSol  = boundary_repair_local(newSol, Xmin, Xmax, 'random');

            remFE = maxFES - FES;
            val_sol = inf(1, 2);
            nEvalScout = min(2, remFE);
            for s = 1:nEvalScout
                val_sol(s) = obj(newSol(s,:));
                FES = FES + 1;
            end

            [~, min_index] = min(val_sol);
            X(ind,:)   = newSol(min_index(1), :);
            val_X(ind) = val_sol(min_index(1));
        end

        % Memorize best
        [currBest, idxBest] = min(val_X);
        if currBest < val_gBest
            val_gBest = currBest;
            gBest = X(idxBest(1), :);
        end

        curve(end+1,1) = val_gBest; %#ok<AGROW>
    end

    out.best_f  = val_gBest;
    out.best_x  = gBest;
    out.curve   = curve;
    out.fe_used = FES;
end

function cfg = fill_defaults(cfg)
    if ~isfield(cfg, 'pop')   || isempty(cfg.pop),   cfg.pop   = 50; end
    if ~isfield(cfg, 'dim')   || isempty(cfg.dim),   cfg.dim   = 20; end
    if ~isfield(cfg, 'maxFE') || isempty(cfg.maxFE), cfg.maxFE = 10000 * cfg.dim; end
    if ~isfield(cfg, 'lb')    || isempty(cfg.lb),    cfg.lb    = -100; end
    if ~isfield(cfg, 'ub')    || isempty(cfg.ub),    cfg.ub    = 100; end
    if ~isfield(cfg, 'limit') || isempty(cfg.limit), cfg.limit = 200; end
    if ~isfield(cfg, 'CR')    || isempty(cfg.CR),    cfg.CR    = 0.5; end
end

function b = expand_bound(v, D)
    if isscalar(v)
        b = ones(1, D) * v;
    else
        b = reshape(v, 1, []);
    end
end

function y = cec22_eval_row(x, fid)
    x = x(:).';
    ok = false;
    try
        y = cec22_test_func(x', fid);
        ok = true;
    catch
    end
    if ~ok
        y = cec22_test_func(x, fid);
    end
    if numel(y) > 1
        y = y(1);
    end
end

function r = generateR_local(popsize, i)
    r1 = randi(popsize);
    while r1 == i
        r1 = randi(popsize);
    end
    r2 = randi(popsize);
    while r2 == r1 || r2 == i
        r2 = randi(popsize);
    end
    r3 = randi(popsize);
    while r3 == r2 || r3 == r1 || r3 == i
        r3 = randi(popsize);
    end
    r4 = randi(popsize);
    while r4 == r3 || r4 == r2 || r4 == r1 || r4 == i
        r4 = randi(popsize);
    end
    r5 = randi(popsize);
    while r5 == r4 || r5 == r3 || r5 == r2 || r5 == r1 || r5 == i
        r5 = randi(popsize);
    end
    r = [r1 r2 r3 r4 r5];
end

function u = boundary_repair_local(v, low, up, str)
    [NP, D] = size(v);
    u = v;

    switch lower(str)
        case 'absorb'
            for i = 1:NP
                for j = 1:D
                    if v(i,j) > up(j)
                        u(i,j) = up(j);
                    elseif v(i,j) < low(j)
                        u(i,j) = low(j);
                    else
                        u(i,j) = v(i,j);
                    end
                end
            end

        case 'random'
            for i = 1:NP
                for j = 1:D
                    if v(i,j) > up(j) || v(i,j) < low(j)
                        u(i,j) = low(j) + rand*(up(j)-low(j));
                    else
                        u(i,j) = v(i,j);
                    end
                end
            end

        case 'reflect'
            for i = 1:NP
                for j = 1:D
                    if v(i,j) > up(j)
                        u(i,j) = max(2*up(j)-v(i,j), low(j));
                    elseif v(i,j) < low(j)
                        u(i,j) = min(2*low(j)-v(i,j), up(j));
                    else
                        u(i,j) = v(i,j);
                    end
                end
            end

        otherwise
            error('Unknown boundary repair strategy: %s', str);
    end
end

function fFitness = calculateFitness_local(fObjV)
    fFitness = zeros(size(fObjV));
    ind = find(fObjV >= 0);
    fFitness(ind) = 1 ./ (fObjV(ind) + 1);
    ind = find(fObjV < 0);
    fFitness(ind) = 1 + abs(fObjV(ind));
end

function index = fitnessDistanceBalance_local(population, fitness)
    [~, bestIndex] = min(fitness);
    best = population(bestIndex, :);
    [populationSize, dimension] = size(population);

    distances = zeros(1, populationSize);
    normFitness = zeros(1, populationSize);
    normDistances = zeros(1, populationSize);

    if min(fitness) == max(fitness)
        index = randi(populationSize);
        return;
    end

    for i = 1:populationSize
        value = 0;
        for j = 1:dimension
            value = value + abs(best(j) - population(i, j));
        end
        distances(i) = value;
    end

    minFitness = min(fitness);
    maxMinFitness = max(fitness) - minFitness;
    minDistance = min(distances);
    maxMinDistance = max(distances) - minDistance;

    if maxMinFitness <= eps
        normFitness(:) = 1;
    else
        for i = 1:populationSize
            normFitness(i) = 1 - ((fitness(i) - minFitness) / maxMinFitness);
        end
    end

    if maxMinDistance <= eps
        normDistances(:) = 0;
    else
        for i = 1:populationSize
            normDistances(i) = (distances(i) - minDistance) / maxMinDistance;
        end
    end

    divDistances = normFitness + normDistances;
    [~, index] = max(divDistances);
end
