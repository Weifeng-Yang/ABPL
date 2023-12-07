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
    [V,L(j),Xtemp,temp]=gradibpg(vv,varz,var,ngmar,j,num,r);
    [var{j},~]=PROX(varze{j},var{j},V,aa(j));
    tao=1/(r*L(j));    
    for k=1:len-1
        vv{j}=var{j}+bt(j)*(var{j}-varK{j});
        varz{j}=var{j}+bt(j)*(var{j}-varK{j});
        varK{j}=var{j};
        mar=varz{j}*temp-Xtemp;
        V=vv{j}-tao*mar;
        [var{j},~]=PROX(varze{j},var{j},V,aa(j));
    end
    end
end
