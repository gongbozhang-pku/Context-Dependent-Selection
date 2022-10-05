function [PCS] = pAOAPmc(k, q, n0, T, mu0, sigma0, sigma, num, m)

PCS1=zeros(T,q,num);

parfor_progress(num);

tic

parfor t=1:num

    truemu = normrnd(mu0, sqrt(sigma0));

    X0 = normrnd(truemu .* ones(k, q, n0), sqrt(sigma) .* ones(k, q, n0));
    estmean = mean(X0, 3);
    estvar = var(X0, 0, 3);

    N=n0 * ones(k, q);

    pv = (1 ./ sigma0 + N ./ estvar) .^ (-1);
    pm = pv .* (mu0 ./ sigma0 + N .* estmean ./ estvar);

    Nrv=normrnd(0, 1, 1, T);

    [~, rbb] = sort(truemu, 1, 'descend');
    rb = sort(rbb(1:m, :), 1);

    for budget=1:T

        [~, id1] = sort(pm, 1, 'descend');
        id11 = sort(id1(1:m, :), 1);

        PCS1(budget, :, t) = prod(id11 == rb) / num;

        col = mod(budget,q);
        col(col == 0) = q;

        cm=setdiff((1:k),id11(:,col));

        V = zeros(1,k);

        for i=1:k
            nv=pv(:,col);
            M=N(:,col);
            M(i)=N(i,col)+1;
            nv(i)=(1/sigma0(i,col)+M(i)/estvar(i,col))^(-1);
            a = 1:m;
            b = 1:(k-m);
            V1 = (pm((id11(a,col)),col)-pm(cm(b),col)').^2./(nv(id11(a,col))+nv(cm(b))');
            V(i) = min(min(V1));
        end
        [~,id2]=max(V);

        mm=estmean(id2,col);
        x=truemu(id2,col)+(sigma(id2,col)).^(1/2).*Nrv(budget);
        estmean(id2,col)=(estmean(id2,col).*N(id2,col)+x)./(N(id2,col)+1);
        estvar(id2,col)=((N(id2,col)-1)*estvar(id2,col)+(x-mm)*(x-estmean(id2,col)))/N(id2,col);

        pv(id2,col)=(1./sigma0(id2,col)+N(id2,col)./estvar(id2,col))^(-1);
        pm(id2,col)=pv(id2,col).*(mu0(id2,col)./sigma0(id2,col)+N(id2,col).*estmean(id2,col)./estvar(id2,col));
    end

    parfor_progress;
end
PCS = min(sum(PCS1, 3), [], 2);
parfor_progress(0);
toc
end