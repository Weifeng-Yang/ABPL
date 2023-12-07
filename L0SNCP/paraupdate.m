function var=paraupdate(var,varze,num,ngmar,aa,r)

for j=1:num
V=grad1(var,ngmar,r,j,num);
[var{j},~]=PROX(varze{j},var{j},V,aa(j));
end


end