%% Gradient calculation function for TiBPALM
function U=gradtibpalm(var,varK,varK1,at,bt,ngmar,n,num,r)
  


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
    ctemp=norm(ZT*ZT','fro');
%     r=2*(at+bt)/ctemp+1;
    ck=2*(at+bt)+ctemp;
    ck=r*checkck(ck);
    mar=var{n}*ZT;
    U=var{n}-1/ck*((mar-ngmar)*ZT'+at*(varK{n}-var{n})+bt*(varK1{n}-varK{n}));  
    elseif(n<num)
    ctemp=norm(Z'*Z,'fro')*norm(ZT*ZT','fro');
    ck=2*(at+bt)+ctemp;
    ck=r*checkck(ck);
    mar=Z*var{n}*ZT;
    U=var{n}-1/ck*(Z'*(mar-ngmar)*(ZT)'+at*(varK{n}-var{n})+bt*(varK1{n}-varK{n}));  
    elseif(n==num)
    ctemp=norm(Z'*Z,'fro');
    ck=2*(at+bt)+ctemp;
    ck=r*checkck(ck);
    mar=Z*var{n};
    U=var{n}-1/ck*(Z'*(mar-ngmar)+at*(varK{n}-var{n})+bt*(varK1{n}-varK{n}));  
    end


end

