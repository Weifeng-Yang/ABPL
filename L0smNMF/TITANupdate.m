%% The parameters update function of TITAN
function [var,L,bt]=TITANupdate(var,varze,varK,num,ngmar,aa,randsor,L,LK,r,wk)
    vv=var;
    bt=zeros(1,num);
    m=unique(randsor,'stable');
    for i=1:length(m)
    len=length(find(randsor==m(i)));
    j=m(i);
    %% This is the extrapolation parameter update strategy adopted by the TITAN algorithm.
    bt(j)=min(wk,(r-1)/(2*r)*sqrt(LK(j)/L(j)));
    %% Update parameters
    vv{j}=var{j}+bt(j)*(var{j}-varK{j});
    varK{j}=var{j};
    [V,L(j),Z,ZT]=gradit(vv,vv,var,ngmar,j,num,r);
    [var{j},~]=PROX(varze{j},V,aa(j));
    tao=1/(r*L(j));    
    for k=1:len-1
        vv{j}=var{j}+bt(j)*(var{j}-varK{j});
        varK{j}=var{j};
        mar=Z*vv{j}*ZT;
        V=vv{j}-tao*Z'*(mar-ngmar)*(ZT)';
        [var{j},~]=PROX(varze{j},V,aa(j));
    end
    end
    
end
