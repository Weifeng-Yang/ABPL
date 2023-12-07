function [var,L]=BPLupdate(var,varze,num,ngmar,aa,r)
randsor=randperm(num);
% L=zeros(num,1);
for i=1:num
% if(randsor(i)==1)
% [U,L(1)]=gradbpl1(var,ngmar,1,num,r);
% [var{1},~]=PROX(varze{1},var{1},U,aa(1));
% elseif(randsor(i)<num)
j=randsor(i);
[V,L(j)]=gradbpl(var,ngmar,j,num,r);
[var{j},~]=PROX(varze{j},var{j},V,aa(j));
% else
% [G,L(num)]=gradbplend(var,ngmar,num,num,r);
% [var{num},~]=PROX(varze{num},var{end},G,aa(end));
% end

end

end