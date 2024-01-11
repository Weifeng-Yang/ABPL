%%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions
%%  This algorithm comes from the paper: "An Inertial Block Majorization Minimization Framework for Nonsmooth Nonconvex Optimization"
function [var,loss,timerun,bts]=TITAN(var,ngmar,aa,maxiteropt,num,stopindex,r)
%% initialization algorithm
loss=[];
timerun=[0];
for i=1:num
    varze{i}=zeros(size(var{i},1),size(var{i},2));
end
bts=[];
LK=zeros(1,num);
L=ones(1,num);
tk=1;
varK=var;
wk=(tk-1)/(tk+1);
loss(1)=compute(var,num,ngmar);
returnloss=norm(ngmar,"fro");


t1=clock;
for i=1:maxiteropt
%% update parameters
fprintf("%d\n",i);
randsor=randperm(num);
randlen=2; %% The number of repeated updates adopts the default value of the TITAN algorithm, which is one time.
randsor=repmat(randsor,1,randlen);
LtempK=LK;
vartempK=varK;
varK=var;
LK=L;
[var,L,bt]=TITANupdate(var,varze,vartempK,num,ngmar,aa,randsor,L,LtempK,r,wk);
loss(i+1)=compute(var,num,ngmar);
bts{i}=bt;


%% Check if termination condition is met
fprintf("TITAN\n");
for j=1:num
    fprintf("nonzero:%d\n",nnz(var{j}));
end
t2=clock;
timerun(i+1)=etime(t2,t1);
ReeK=loss(i)/returnloss;
Ree=loss(i+1)/returnloss;
Res=abs(Ree-ReeK);
fprintf("Rel: %d\n",Res);
stop=stopcheck(timerun,stopindex);
if(stop==1)
    fprintf("Number of terminations：%d\n",i);
    pause(4);
    break;
end
tk=(1+sqrt(1+4*tk^2))/2;
wk=(tk-1)/(tk+1);
end
end













