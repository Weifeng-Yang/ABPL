%% Gradient calculation function for IBPG and TITAN
function [U,L,Z,ZT]=gradit(var,varz,varK,ngmar,n,num,r)
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
    ck=checkck(ck);
    mar=varz{n}*ZT;
    U=var{n}-1/(r*ck)*(mar-ngmar)*(ZT)';   
    L=ck;
    elseif(n<num)
    ck=norm(Z'*Z,'fro')*norm(ZT*ZT','fro');
    ck=checkck(ck);
    L=ck;
    mar=Z*varz{n}*ZT;
    U=var{n}-1/(r*ck)*Z'*(mar-ngmar)*(ZT)';   
    elseif(n==num)
    ck=norm(Z'*Z,'fro');
    ck=checkck(ck);
    L=ck;
    ZT=eye(size(var{n},2),size(var{n},2)); 
    mar=Z*varz{n};
    U=var{n}-1/(r*ck)*Z'*(mar-ngmar);   
    end 
end
