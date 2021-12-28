function [PCS,N]=AOAPmc(k,q,n0,T,mu0,sigma0,sigma,num,m)

PCS1=zeros(T,q);
estmean = zeros(k,q);
estvar = zeros(k,q);

mu=zeros(k,q);

rb=zeros(k,q);

X0 = zeros(n0,q);

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
        estvar(i,:)=var(X0); 
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
        
        cm=setdiff((1:k),id11(:,col));
        
        for i=1:k  
            nv=pv(:,col);
            M=N(:,col);  
            M(i)=N(i,col)+1; 
            nv(i)=(1/sigma0(i,col)+M(i)/estvar(i,col))^(-1); 
            for a=1:m
                for b=1:(k-m)
                    V1(a,b)= (pm((id11(a,col)),col)-pm(cm(b),col)).^2/(nv(id11(a,col))+nv(cm(b)));
                end
            end
            V(i)=min(min(V1));
        end
        [~,id2]=max(V);
        
        mm=estmean(id2,col);
        x=mu(id2,col)+(sigma(id2,col)).^(1/2).*Nrv(budget); 
        estmean(id2,col)=(estmean(id2,col).*N(id2,col)+x)./(N(id2,col)+1);
        estvar(id2,col)=(N(id2,col)./(N(id2,col)+1)).*(estvar(id2,col)+(mm-x).^2./(N(id2,col)+1)); 
        N(id2,col)=N(id2,col)+1;
        
        pv(id2,col)=(1./sigma0(id2,col)+N(id2,col)./estvar(id2,col))^(-1);
        pm(id2,col)=pv(id2,col).*(mu0(id2,col)./sigma0(id2,col)+N(id2,col).*estmean(id2,col)./estvar(id2,col)); 
    end
end
toc
PCS = min((PCS1)');
end
