function var=paraupdate(var,varze,num,ngmar,aa,r)



for j=1:num

V{j}=grad(var,ngmar,r,j,num);
[var{j},~]=PROX(varze{j},var{j},V{j},aa(j));

end




end