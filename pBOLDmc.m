function [PCS] = pBOLDmc(k, q, n0, T, mu0, sigma0, sigma, num, m)

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

    wt = ones(k,q) / (k*q);

    for budget = 1:T

        [~, id1] = sort(pm, 1, 'descend');
        id11 = sort(id1(1:m, :), 1);

        PCS1(budget, :, t) = prod(id11 == rb) / num;

        [~, id3] = sort(estmean, 1, 'descend');
        id33 = sort(id3(1:m, :), 1);
        cm = sort(id3(m+1:k, :), 1);

        V = zeros(3,q);
        for jj = 1:q
            a = 1:m;
            b = 1:(k-m);
            V1 = (estmean(id33(a,jj),jj)-estmean(cm(b,jj),jj)').^2./(estvar(id33(a,jj),jj)./wt(id33(a,jj),jj)+estvar(cm(b,jj),jj)'./wt(cm(b,jj),jj)');
            V(1,jj) = min(min(V1));
            [A,B] = find(V1 == min(min(V1)));
            V(2,jj) = A;
            V(3,jj) = B;
        end
        [~,id22] = min(V(1,:));

        u0 = sum(N(id33(:,id22),id22).^2./estvar(id33(:,id22),id22));
        un = sum(N(cm(:,id22),id22).^2./estvar(cm(:,id22),id22));
        if u0<un
            id21=id33(V(2,id22),id22);
        else
            id21=cm(V(3,id22),id22);
        end

        mm = estmean(id21, id22);
        x = truemu(id21, id22) + (sigma(id21, id22)) .^ (1/2) .* Nrv(budget);
        estmean(id21, id22) = (estmean(id21, id22) .* N(id21, id22) + x) ./ (N(id21, id22) + 1);
        estvar(id21, id22)=((N(id21, id22)-1)*estvar(id21, id22)+(x-mm)*(x-estmean(id21, id22)))/N(id21, id22);
        N(id21, id22) = N(id21, id22) + 1;

        wt = N./sum(sum(N));

        pv(id21, id22) = (1 ./ sigma0(id21, id22) + N(id21, id22) ./ estvar(id21, id22)) .^ (-1);
        pm(id21, id22) = pv(id21, id22) .* (mu0(id21, id22) ./ sigma0(id21, id22) + N(id21, id22) .* estmean(id21, id22) ./ estvar(id21, id22));
    end

    parfor_progress;
end
PCS = min(sum(PCS1, 3), [], 2);
parfor_progress(0);
toc
end