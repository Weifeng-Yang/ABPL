%% Check whether the decomposed matrix meets the L0 norm requirement
function flag=chencknon(w,aa,num)
flag=0;
 for i=1:num
     ab=nnz(w{i}~=0);
     as=nnz(w{i}<0);
     if(as>0 || ab>aa(i))
         flag=1;
         return;
     end
 end
end