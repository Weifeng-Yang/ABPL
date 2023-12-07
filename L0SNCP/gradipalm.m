%% Gradient calculation function for algorithms iPALM
function U=gradipalm(var,varz,varK,at,bt,ngmar,n,num)
%    B=krob1(varK,n,num);
  [Xtemp,temp]=krob2(varK,n,num,ngmar);

%    temp=B'*B;
%    Xtemp=X*B;
   a=0;
   ck=(1+2*bt)/(1-2*at)*norm(temp,'fro');
   tao=1/ck;
   if(temp==0) 
       while(a==0)
           a=rand(1);
       end
       tao=1/a;
   end

   mar=varz{n}*temp-Xtemp;
   U=var{n}-tao*mar;


end