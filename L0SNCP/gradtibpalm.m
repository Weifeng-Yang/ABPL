%% Gradient calculation function for iPALM
function U=gradtibpalm(var,varK,varK1,at,bt,ngmar,n,num,r)
  [Xtemp,temp]=krob2(varK,n,num,ngmar);
   a=0;
   ctemp=norm(temp,'fro');
   ck=2*(at+bt)+ctemp;
   tao=1/(r*ck);
   mar=var{n}*temp-Xtemp;
   U=var{n}-tao*(mar+at*(varK{n}-var{n})+bt*(varK1{n}-varK{n}));
end