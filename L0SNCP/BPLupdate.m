function [var,L]=BPLupdate(var,varze,num,ngmar,aa,r)


randsor=randperm(num);
for i=1:num
j=randsor(i);
[V,L(j)]=gradibpl(var,ngmar,j,num,r);
[var{j},~]=PROX(varze{j},var{j},V,aa(j));
end


end