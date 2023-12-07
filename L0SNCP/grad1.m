%% Gradient calculation function for algorithms PALM, APGnc, ABPL
function U=grad1(var,ngmar,r,n,num)
%    B=krob1(var,n,num);
  [Xtemp,temp]=krob2(var,n,num,ngmar);

%    temp1=B'*B;
%    Xtemp1=X*B;

   a=0;
   ck=r*norm(temp,'fro');
   tao=1/ck;
   if(temp==0)
       while(a==0)
           a=rand(1);
       end
       tao=1/a;
   end

   mar=var{n}*temp-Xtemp;
   U=var{n}-tao*mar;


end