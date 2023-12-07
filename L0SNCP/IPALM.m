%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions

function [var,loss,timerun,bts,ats,rate]=IPALM(var,ngmar,aa,maxiteropt,bt,obj,num,N,stopindex,r)
%% initialization algorithm
loss=[];
timerun=[0];
for i=1:num
    varze{i}=zeros(size(ngmar,i),N);
end
cardinal=1e-8;
bts=[];
rate=[0];
varK=var;
for j=1:num
    fprintf("nonzero:%d\n",nnz(var{j}~=0));
end
returnloss=norm(ngmar);
loss(1)=compute(var,ngmar);
t1=clock;


for i=1:maxiteropt
%% update parameters
at=(i-1)/(2*(i+2));
bt=(i-1)/(i+2);  
fprintf("%d\n",i);

for j=1:num   
    vv{j}=var{j}+at*(var{j}-varK{j});
    varz{j}=var{j}+bt*(var{j}-varK{j});
end
varK=var;
var=IPALMupdate(vv,varze,varK,varz,num,ngmar,aa,at,bt);
loss(i+1)=compute(var,ngmar);
bts(i)=bt;
ats(i)=at;



%% Check if termination condition is met
fprintf("iPALM\n");
for j=1:num
    fprintf("nonzero:%d\n",nnz(var{j}~=0));
%     check=check+norm(var{j}-varK{j},'fro');    
end

t2=clock;
timerun(i+1)=etime(t2,t1);
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














