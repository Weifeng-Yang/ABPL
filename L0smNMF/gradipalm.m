%% Gradient calculation function for iPALM
function U=gradipalm(var,varz,varK,at,bt,ngmar,n,num)
    if(n>1)
    Z=eye(size(var{1},1),size(var{1},1));
    end
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
    ck=(1+2*bt)/(1-2*at)*norm(ZT*ZT','fro');
    ck=checkck(ck);
    mar=varz{n}*ZT;
    U=var{n}-1/ck*(mar-ngmar)*ZT';  
    elseif(n<num)
    ck=(1+2*bt)/(1-2*at)*norm(Z'*Z,'fro')*norm(ZT*ZT','fro');
    ck=checkck(ck);
    mar=Z*varz{n}*ZT;
    U=var{n}-1/ck*Z'*(mar-ngmar)*(ZT)';  
    elseif(n==num)
    ck=(1+2*bt)/(1-2*at)*norm(Z'*Z,'fro');
    ck=checkck(ck);
    mar=Z*varz{n};
    U=var{n}-1/ck*Z'*(mar-ngmar);  
    end


end

