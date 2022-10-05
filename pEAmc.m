function [PCS] = pEAmc(k, q, n0, T, mu0, sigma0, sigma, num, m)

PCS1 = zeros(T, q, num);

parfor_progress(num);

tic

parfor t = 1:num

    truemu = normrnd(mu0, sqrt(sigma0));

    X0 = normrnd(truemu .* ones(k, q, n0), sqrt(sigma) .* ones(k, q, n0));
    estmean = mean(X0, 3);
    estvar = var(X0, 0, 3);

    N = n0 * ones(k, q);

    pv = (1 ./ sigma0 + N ./ estvar) .^ (-1);
    pm = pv .* (mu0 ./ sigma0 + N .* estmean ./ estvar);

    Nrv = normrnd(0, 1, 1, T);

    [~, rbb] = sort(truemu, 1, 'descend');
    rb = sort(rbb(1:m, :), 1);

    for budget = 1:T

        [~, id1] = sort(pm, 1, 'descend');
        id11 = sort(id1(1:m, :), 1);

        PCS1(budget, :, t) = prod(id11 == rb) / num;

        y1 = mod(budget, k * q);
        if y1 == 0
            id22 = q;
            id21 = k;
        elseif mod(y1, q) == 0
            id22 = q;
            id21 = fix(y1 / q);
        else
            id22 = mod(y1, q);
            id21 = ceil(y1 / q);
        end

        mm = estmean(id21, id22);
        x = truemu(id21, id22) + (sigma(id21, id22)) .^ (1/2) .* Nrv(budget);
        estmean(id21, id22) = (estmean(id21, id22) .* N(id21, id22) + x) ./ (N(id21, id22) + 1);
        estvar(id21, id22)=((N(id21, id22)-1)*estvar(id21, id22)+(x-mm)*(x-estmean(id21, id22)))/N(id21, id22);
        N(id21, id22) = N(id21, id22) + 1;

        pv(id21, id22) = (1 ./ sigma0(id21, id22) + N(id21, id22) ./ estvar(id21, id22)) .^ (-1);
        pm(id21, id22) = pv(id21, id22) .* (mu0(id21, id22) ./ sigma0(id21, id22) + N(id21, id22) .* estmean(id21, id22) ./ estvar(id21, id22));
    end

    parfor_progress;
end
PCS = min(sum(PCS1, 3), [], 2);
parfor_progress(0);
toc
end