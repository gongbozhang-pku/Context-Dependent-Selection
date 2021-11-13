function [PCS,N]=OCBAmc(k,q,n0,T,mu0,sigma0,sigma,num,m)

PCS1=zeros(T,q);
estmean = zeros(k,q);
estvar = zeros(k,q);

mu=zeros(k,q);

rb=zeros(k,q);

X0 = zeros(n0,q);

bar=waitbar(0,'读取数据中...');
tic
for t=1:num
    
    for i=1:k
        mu(i,:)=normrnd(mu0(i,:),(sigma0(i,:)).^(1/2));
    end
    
    for i=1:k
        for j=1:n0
            X0(j,:)=normrnd(mu(i,:),sigma(i,:).^(1/2));
        end
        estmean(i,:)=mean(X0);
        estvar(i,:)=var(X0,1);
    end
    
    N=n0*ones(k,q);
    
    pv = (1./sigma0+N./estvar).^(-1);
    pm = pv.*(mu0./sigma0+N.*estmean./estvar);
    
    Nrv=normrnd(0,1,1,T);
    
    for i=1:q
        [~,rbb]=sort(mu(:,i),'descend');
        rbb=sort(rbb(1:m));
        %rbb=sort(rbb(1:m(i)));
        rb((1:m),i) = rbb;
        %rb((1:m(i)),i) = rbb;
    end
    
    wt=ones(k,q)/k;
    
    for budget=1:T
        
        id11=zeros(k,q);
        
        for i=1:q
            [~,id1]=sort(pm(:,i),'descend');
            id1 = sort(id1(1:m));
            %id1 = sort(id1(1:m(i)));
            id11((1:m),i) = id1;
            %id11((1:m(i)),i) = id1;
        end
        
        for i=1:q
            if id11(:,i)==rb(:,i)
                PCS1(budget,i)= PCS1(budget,i) + 1/num;
            end
        end
        
        res = mod(budget,q);
        if res == 0
            col = q;
        else
            col = res;
        end
        
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
        x=mu(id2,col)+(sigma(id2,col)).^(1/2).*Nrv(budget);
        estmean(id2,col)=(estmean(id2,col).*N(id2,col)+x)./(N(id2,col)+1);
        estvar(id2,col)=(N(id2,col)./(N(id2,col)+1)).*(estvar(id2,col)+(mm-x).^2./(N(id2,col)+1));
        N(id2,col)=N(id2,col)+1;
        
        wt(:,col)=N(:,col)./sum(N(:,col));
        
        pv(id2,col)=(1./sigma0(id2,col)+N(id2,col)./estvar(id2,col))^(-1);
        pm(id2,col)=pv(id2,col).*(mu0(id2,col)./sigma0(id2,col)+N(id2,col).*estmean(id2,col)./estvar(id2,col));
    end
    str=['计算中...',num2str(100*t/num),'%'];
    waitbar(t/num,bar,str)
end
toc
PCS = min((PCS1)');
close(bar)
end