function [PCS] = pmLinGapE(k, q, n0, T, mu0, sigma0, sigma, num, m, cv0, cvsigma0)

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

    cv = normrnd(cv0, sqrt(cvsigma0));

    for budget = 1:T

        [~, id1] = sort(pm, 1, 'descend');
        id11 = sort(id1(1:m, :), 1);
        cm = sort(id1(m+1:k, :), 1);

        PCS1(budget, :, t) = prod(id11 == rb) / num;


        Vmatrix = sqrt(cvsigma0(1)) / 20 + sum(sum(N).*cv.*cv);

        Theta = sum(estmean.*cv,2)./Vmatrix;

        %         id22 = mod(budget,q);
        %         id22(id22 == 0) = q;


        %         a = 1:m;
        %         b = 1:(k-m);
        %
        %         V1 = pm(id11(a,id22),id22)-pm(cm(b,id22),id22)';
        %
        %         [Bid1,Bid2] = find(V1 == min(min(V1)));
        %
        %         if pv(cm(Bid2,id22),id22) <= pv(id11(Bid1,id22),id22)
        %             id21 = id11(Bid1,id22);
        %         else
        %             id21 = cm(Bid2,id22);
        %         end


        V = zeros(2,q);
        for jj = 1:q
            V1 = (Theta(id11(a,jj))-Theta(cm(b,jj))').*cv(jj);
            [IndV1,~] = find(V1 == max(max(V1)));
            V(1,jj) = max(max(V1));
            V(2,jj) = IndV1;
        end
        [~,Bid2] = min(V(1,:));
        Bid1 = V(2,Bid2);

        V2 = zeros(k-m,q);
        for jj = 1:q
            V2 (:,jj) = Theta(id11(Bid1,Bid2))*cv(Bid2)-Theta(cm(b,jj)).*cv(jj)+sqrt(2*log((log(budget)+1)/0.05))*sqrt((cv(Bid2)-cv(jj))^2*cvsigma0(1)/Vmatrix);
        end
        [Cid1,Cid2] = find(V2 == min(max(V2)));

        if Bid2 == Cid2
            if rand(1) <= 0.5
                id21 = id11(Bid1,Bid2);
                id22 = Bid2;
            else
                id21 = cm(Cid1,Cid2);
                id22 = Cid2;
            end
        else
            if cv(Cid2)^2 <= cv(Bid2)^2
                id21 = id11(Bid1,Bid2);
                id22 = Bid2;
            else
                id21 = cm(Cid1,Cid2);
                id22 = Cid2;
            end
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