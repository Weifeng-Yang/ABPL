function var=paraupdaterd(var,varze,num,ngmar,aa,r)
randsor=randperm(num);
for i=1:num
j=randsor(i);
V=grad1(var,ngmar,r,j,num);
[var{j},~]=PROX(varze{j},var{j},V,aa(j));
end


end