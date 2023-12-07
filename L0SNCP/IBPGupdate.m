function [var,L,at,bt]=IBPGupdate(var,varze,varK,varz,num,ngmar,aa,randsor,L,LK,r,wk)
    vv=var;
    at=zeros(1,num);
    bt=at;
    m=unique(randsor,'stable');
    for i=1:length(m)
    len=length(find(randsor==m(i)));
    j=m(i);
    at(j)=min(wk,0.2*(r-1)/r*sqrt(LK(j)/L(j)));
    bt(j)=min(wk,0.2*(r-1)*sqrt(LK(j)/L(j)));
    vv{j}=var{j}+at(j)*(var{j}-varK{j});
    varz{j}=var{j}+bt(j)*(var{j}-varK{j});
    varK{j}=var{j};
    [V,L(j),Xtemp,temp]=gradibpg(vv,varz,var,ngmar,j,num,r);
    [var{j},~]=PROX(varze{j},var{j},V,aa(j));
    tao=1/(r*L(j));    
    for k=1:len-1
        vv{j}=var{j}+at(j)*(var{j}-varK{j});
        varz{j}=var{j}+bt(j)*(var{j}-varK{j});
        varK{j}=var{j};
        mar=varz{j}*temp-Xtemp;%%注意ipalm中第i项是z_{i}^k!
        V=vv{j}-tao*mar;
        [var{j},~]=PROX(varze{j},var{j},V,aa(j));
    end
    end

end
