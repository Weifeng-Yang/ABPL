%% Check whether the decomposed matrix meets the L0 norm requirement
function flag=chencknon(w,aa,num)
flag=0;
 for i=1:num
     if(nnz(w{i}~=0)>aa(i) || nnz(w{i}<0)>0)
         flag=1;
         return;
     end
 end
end