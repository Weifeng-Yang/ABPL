%% Gradient calculation function for algorithms IBPG and TITAN
function [U,L,Z,ZT]=gradibpg(var,varz,varK,ngmar,n,num,r)
    Z=eye(size(var{1},1),size(var{1},1));
    if(n<num)
    ZT=eye(size(var{n+1},1),size(var{n+1},1)); 
    end

    for i=1:n-1
        Z=Z*varK{i};
    end
    for i=n+1:num
        ZT=ZT*varK{i};
    end
    
    if(n==1)
    ck=norm(ZT*ZT','fro');
    elseif(n<num)
    ck=norm(Z'*Z,'fro')*norm(ZT*ZT','fro');
    elseif(n==num)
    ck=norm(Z'*Z,'fro');
    ZT=eye(size(var{n},2),size(var{n},2)); 
    end
    L=ck;
    tao=1/(r*ck);
    if(ck==0)
        a=0;
       while(a==0)
           a=rand(1);
       end
       ck=a;
     end
    mar=Z*varz{n}*ZT;
    U=var{n}-tao*Z'*(mar-ngmar)*(ZT)';   


end
