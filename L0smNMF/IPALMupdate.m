function var=IPALMupdate(var,varze,varK,varz,num,ngmar,aa,at,bt)



for j=1:num

V=gradipalm(var,varz,varK,at,bt,ngmar,j,num);
[var{j},~]=PROX(varze{j},var{j},V,aa(j));

end




end