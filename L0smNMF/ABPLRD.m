%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions
function [var,loss,timerun,bts,rate]=ABPLRD(var,ngmar,aa,maxiteropt,bt,obj,num,stopindex,r,t,btmax)
%% initialization algorithm
loss=[];
timerun=[0];
tk=t;
tks=1;

rate=[0];
for i=1:num
    varze{i}=zeros(size(var{i},1),size(var{i},2));
end
if(obj==2)
    tktemp=tks;
    tks=(1+sqrt(1+4*tks^2))/2;
    bt=(tktemp-1)/(tks);
end
bts=[];
varK=var;

loss(1)=compute(var,num,ngmar);
cardinal=1e-6;
returnloss=norm(ngmar,"fro");


t1=clock;
for i=1:maxiteropt
%% update parameters
fprintf("%d\n",i);
for j=1:num
    vv{j}=var{j}+bt*(var{j}-varK{j});
end
varK=var;
var=paraupdaterd(vv,varze,num,ngmar,aa,r);
loss(i+1)=compute(var,num,ngmar);

%% Judging whether to extrapolate
[bt,var,loss(i+1)]=ABPLRDjudge(var,num,loss(i+1),loss(i),ngmar,varK,varze,aa,vv,bt,tk,obj,btmax,r);
if(obj==2)
    tktemp=tks;
    tks=(1+sqrt(1+4*tks^2))/2;
    bt=(tktemp-1)/(tks);
end
bts(i)=bt;


%% Check if termination condition is met
fprintf("ABPL-random\n");
t2=clock;
timerun(i+1)=etime(t2,t1);
for j=1:num
    fprintf("nonzero:%d\n",nnz(var{j}~=0));
end
ReeK=loss(i)/returnloss;

Ree=loss(i+1)/returnloss;

Res=abs(Ree-ReeK);
rate(i+1)=Res;
fprintf("cri：%d\n",Res);
stop=stopcheck(timerun,stopindex);
if(stop==1)
    fprintf("Number of terminations：%d\n",i);
    siter=i;
    pause(4);
    break;
end




end

end




