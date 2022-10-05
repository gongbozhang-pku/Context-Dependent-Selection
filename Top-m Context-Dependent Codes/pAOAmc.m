function [PCS]=pAOAmc(k,q,n0,T,mu0,sigma0,sigma,num,m)

PCS1 = zeros(T,q,num);

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
        cm = sort(id1(m+1:k, :), 1);
        
        PCS1(budget, :, t) = prod(id11 == rb) / num;
        
        V = zeros(k,q);
        V2 = zeros(1,q);
        
        for i = 1:k
            for j = 1:q
                nv = pv;
                M = N;
                M(i,j) = N(i,j)+1;
                nv(i,j) = (1/sigma0(i,j)+M(i,j)/estvar(i,j))^(-1);
                for jj = 1:q
                    a = 1:m;
                    b = 1:(k-m);
                    V1 = (pm(id11(a,jj),jj)-pm(cm(b,jj),jj)').^2./(nv(id11(a,jj),jj)+nv(cm(b,jj),jj)');
                    V2(jj) = min(min(V1));
                end
                V(i,j) = min(V2);
            end
        end
        
        [id21,id22] = find(V == max(max(V)));
        
        mm=estmean(id21,id22);
        x=truemu(id21,id22)+(sigma(id21,id22)).^(1/2).*Nrv(budget);
        estmean(id21,id22)=(estmean(id21,id22).*N(id21,id22)+x)./(N(id21,id22)+1);
        estvar(id21,id22)= ((N(id21,id22)-1)*estvar(id21,id22)+(x-mm)*(x-estmean(id21,id22)))/N(id21,id22);
        N(id21,id22)=N(id21,id22)+1;
        
        pv(id21,id22)=(1./sigma0(id21,id22)+N(id21,id22)./estvar(id21,id22)).^(-1);
        pm(id21,id22)=pv(id21,id22).*(mu0(id21,id22)./sigma0(id21,id22)+N(id21,id22).*estmean(id21,id22)./estvar(id21,id22));
    end
    
    parfor_progress;
end
PCS = min(sum(PCS1, 3), [], 2);
parfor_progress(0);
toc
end