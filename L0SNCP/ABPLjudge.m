function [bt,var,loss]=ABPLjudge(var,num,nga,ngak,ngmar,aa,vv,bt,tk,obj,varK,varze,r,btmax)
flag=0;
index=chencknon(vv,aa,num);
checkside=0;

if(index==0)
ngaV=compute(vv,ngmar);
    if(ngaV<ngak)
        flag=1;
    else
        checkside=1;
    end
else
    checkside=1;
end

if(checkside)
    vartemp=paraupdate(varK,varze,num,ngmar,aa,r);
    ngai=compute(vartemp,ngmar);
    if(nga<ngai)
        flag=1;
    else
        flag=0;
   end

end

if(flag==1)         
        if(obj==1)
            bt=min(bt*tk,btmax);
        end
        loss=nga;
else    
        var=vartemp;
        if(obj==1)
            bt=bt/tk; 
        end
        loss=ngai;

end



end






