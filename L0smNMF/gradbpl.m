%% Gradient calculation function for algorithms BPL
function [U,L]=gradbpl(var,ngmar,n,num,r)
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
    ck=norm(ZT*ZT','fro');
    L=ck;
    ck=checkck(ck);
    mar=var{n}*ZT;
    U=var{n}-1/(r*ck)*(mar-ngmar)*ZT'; 
    elseif(n<num)
    ck=norm(Z'*Z,'fro')*norm(ZT*ZT','fro');
    L=ck;
    ck=checkck(ck);
    mar=Z*var{n}*ZT;
    U=var{n}-1/(r*ck)*Z'*(mar-ngmar)*(ZT)';  
    elseif(n==num)
    ck=norm(Z'*Z,'fro');
    L=ck;
    ck=checkck(ck);
%     ZT=eye(size(var{n},2),size(var{n},2)); 
    mar=Z*var{n};
    U=var{n}-1/(r*ck)*Z'*(mar-ngmar);   
    end


end


