%% Gradient calculation function for algorithms BPL
function [U,L]=gradbpl(var,ngmar,n,num,r)
%    B=krob1(var,n,num);
  [Xtemp,temp]=krob2(var,n,num,ngmar);

%    temp=B'*B;
%    Xtemp=X*B;

   a=0;
   ck=norm(temp,'fro');
   L=ck;
   tao=1/(r*ck);
   if(temp==0) 
       while(a==0)
           a=rand(1);
       end
       tao=1/a;
   end

   mar=var{n}*temp-Xtemp;
   U=var{n}-tao*mar;


end