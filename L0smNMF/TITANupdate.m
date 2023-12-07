function [var,L,bt]=TITANupdate(var,varze,varK,varz,num,ngmar,aa,randsor,L,LK,r,wk)
   vv=var;
    bt=zeros(1,num);
    m=unique(randsor,'stable');
    
    for i=1:length(m)
    len=length(find(randsor==m(i)));
    j=m(i);
    bt(j)=min(wk,(r-1)/(2*r)*sqrt(LK(j)/L(j)));
    vv{j}=var{j}+bt(j)*(var{j}-varK{j});
    varz{j}=var{j}+bt(j)*(var{j}-varK{j});
    varK{j}=var{j};
    [V,L(j),Z,ZT]=gradibpg(vv,varz,var,ngmar,j,num,r);
    [var{j},~]=PROX(varze{j},var{j},V,aa(j));
    tao=1/(r*L(j));    
    for k=1:len-1
        vv{j}=var{j}+bt(j)*(var{j}-varK{j});
        varz{j}=var{j}+bt(j)*(var{j}-varK{j});
        varK{j}=var{j};
        mar=Z*varz{j}*ZT;
        V=vv{j}-tao*Z'*(mar-ngmar)*(ZT)';

        [var{j},~]=PROX(varze{j},var{j},V,aa(j));
    end
    end


% G=gradipalmend(var,varz,varK,at,bt,ngmar,num,num);
% [var{num},~]=PROX(varze{num},var{end},G,aa(end));
% varK{num}=var{num};

end
