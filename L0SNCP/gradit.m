%% Gradient calculation function for IBPG and TITAN
function [U,L,Xtemp,temp]=gradit(var,varz,varK,ngmar,n,num,r)
  [Xtemp,temp]=krob2(varK,n,num,ngmar);
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
   mar=varz{n}*temp-Xtemp;
   U=var{n}-tao*mar;
end