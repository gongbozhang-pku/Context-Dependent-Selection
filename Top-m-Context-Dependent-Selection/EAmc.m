function [PCS,N]=EAmc(k,q,n0,T,mu0,sigma0,sigma,num,m)

tic
PCS1 = zeros(T,q);
estmean = zeros(k,q);
estvar = zeros(k,q);

mu=zeros(k,q);
rb=zeros(m,q);
X0 = zeros(n0,q);

bar=waitbar(0,'读取数据中...');

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
        rb(:,i) = rbb;
    end
    
    for budget=1:T
        
        id11=zeros(m,q);
        cm=zeros((k-m),q);
        for i=1:q
            [~,id1]=sort(pm(:,i),'descend');
            id1 = sort(id1(1:m));
            id11(:,i) = id1;
        end
        
        for i=1:q      
            if id11(:,i)==rb(:,i)
                PCS1(budget,i)= PCS1(budget,i) + 1/num;
            end
        end
        
        y1 = mod(budget,k*q);
        if y1 == 0
            id22 = q;
            id21 = k;
        elseif mod(y1,q)==0
            id22 = q;
            id21 = fix(y1/q);
        else
            id22 = mod(y1,q);
            id21 = ceil(y1/q);
        end
        
        mm=estmean(id21,id22);
        x=mu(id21,id22)+(sigma(id21,id22)).^(1/2).*Nrv(budget);
        estmean(id21,id22)=(estmean(id21,id22).*N(id21,id22)+x)./(N(id21,id22)+1);
        estvar(id21,id22)=(N(id21,id22)./(N(id21,id22)+1)).*(estvar(id21,id22)+(mm-x).^2./(N(id21,id22)+1));
        N(id21,id22)=N(id21,id22)+1;
        
        pv(id21,id22)=(1./sigma0(id21,id22)+N(id21,id22)./estvar(id21,id22)).^(-1);
        pm(id21,id22)=pv(id21,id22).*(mu0(id21,id22)./sigma0(id21,id22)+N(id21,id22).*estmean(id21,id22)./estvar(id21,id22));
    end
    str=['计算中...',num2str(100*t/num),'%'];
    waitbar(t/num,bar,str)
end
PCS = min((PCS1)');
close(bar)
toc 
end