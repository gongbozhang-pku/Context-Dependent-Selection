function [PCS]=pOCBAmc(k,q,n0,T,mu0,sigma0,sigma,num,m)

PCS1=zeros(T,q,num);

parfor_progress(num);

tic

parfor t=1:num

    truemu = normrnd(mu0, sqrt(sigma0));
    
    X0 = normrnd(truemu .* ones(k, q, n0), sqrt(sigma) .* ones(k, q, n0));
    estmean = mean(X0, 3);
    estvar = var(X0, 0, 3);
    
    N = n0 * ones(k,q);
    
    pv = (1 ./ sigma0 + N ./ estvar) .^ (-1);
    pm = pv .* (mu0 ./ sigma0 + N .* estmean ./ estvar);
    
    Nrv = normrnd(0, 1, 1, T);
    
    [~, rbb] = sort(truemu, 1, 'descend');
    rb = sort(rbb(1:m, :), 1);
    
    wt = ones(k,q) / k;

    for budget=1:T
        
        [~, id1] = sort(pm, 1, 'descend');
        id11 = sort(id1(1:m, :), 1);
        
        PCS1(budget, :, t) = prod(id11 == rb) / num;
        
        col = mod(budget,q);
        col(col == 0) = q;
        
        [pmm,id4]=sort(estmean(:,col),'descend');
        c=(estvar(id4(m+1),col)*pmm(m)+estvar(id4(m),col)*pmm(m+1))/(estvar(id4(m),col)+estvar(id4(m+1),col));
        
        w = zeros(1,k);
        delta = estmean(:,col)-c;
        Omega=1:k;
        AW1=(estvar(Omega,col).*delta(k).^2)./(estvar(k,col).*delta(Omega).^2);
        w(k)=1/sum(AW1);
        w(Omega)=AW1*w(k);
        df=(k*n0+budget)*w'-(k*n0+budget-1)*wt(:,col);
        
        [~,id2]=max(df);
        
        mm=estmean(id2,col);
        x=truemu(id2,col)+(sigma(id2,col)).^(1/2).*Nrv(budget);
        estmean(id2,col)=(estmean(id2,col).*N(id2,col)+x)./(N(id2,col)+1);
        estvar(id2,col)=((N(id2,col)-1)*estvar(id2,col)+(x-mm)*(x-estmean(id2,col)))/N(id2,col);
        N(id2,col)=N(id2,col)+1;

        wt(:,col)=N(:,col)./sum(N(:,col));
        
        pv(id2,col)=(1./sigma0(id2,col)+N(id2,col)./estvar(id2,col))^(-1);
        pm(id2,col)=pv(id2,col).*(mu0(id2,col)./sigma0(id2,col)+N(id2,col).*estmean(id2,col)./estvar(id2,col));
    end
    
    parfor_progress;
end
PCS = min(sum(PCS1, 3), [], 2);
parfor_progress(0);
toc
end