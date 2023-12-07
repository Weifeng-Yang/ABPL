%% Gradient calculation function for algorithms PALM, APGnc, ABPL
function U=grad(var,ngmar,r,n,num)
    if(n>1)
    Z=eye(size(var{1},1),size(var{1},1));
    end
    if(n<num)
    ZT=eye(size(var{n+1},1),size(var{n+1},1)); 
    end

    for i=1:n-1
        Z=Z*var{i};
    end
    for i=n+1:num
        ZT=ZT*var{i};
    end

    if(n==1)
    ck=r*norm(ZT*ZT','fro');
    ck=checkck(ck);
    mar=var{n}*ZT;
    U=var{n}-1/ck*(mar-ngmar)*ZT';  
    elseif(n<num)
    ck=r*norm(Z'*Z,'fro')*norm(ZT*ZT','fro');
    ck=checkck(ck);
    mar=Z*var{n}*ZT;
    U=var{n}-1/ck*Z'*(mar-ngmar)*(ZT)';  
    elseif(n==num)
    ck=r*norm(Z'*Z,'fro');
%     ZT=eye(size(var{n},2),size(var{n},2)); 
    ck=checkck(ck);
    mar=Z*var{n};
    U=var{n}-1/ck*Z'*(mar-ngmar);   
    end


end

