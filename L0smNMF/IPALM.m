%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions

function [var,loss,timerun,bts,ats,rate]=IPALM(var,ngmar,aa,maxiteropt,bt,num,stopindex)
%% initialization algorithm
loss=[];
timerun=[0];
rate=[0];
for i=1:num
    varze{i}=zeros(size(var{i},1),size(var{i},2));
end
bts=[];
varK=var;


loss(1)=compute(var,num,ngmar);

cardinal=1e-6;
returnloss=norm(ngmar,"fro");

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
loss(i+1)=compute(var,num,ngmar);
bts(i)=bt;
ats(i)=at;





%% Check if termination condition is met
fprintf("iPALM\n");
for j=1:num
    fprintf("nonzero:%d\n",nnz(var{j}));  
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












function X=Xcomp(ngmar,num)
ngmar=double(ngmar);
X={};
for i=1:num
X{i}=ten2mat(ngmar,size(ngmar),i);
end
end


