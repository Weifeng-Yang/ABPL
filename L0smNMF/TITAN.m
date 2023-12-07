%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions

function [var,loss,timerun,bts,rate]=TITAN(var,ngmar,aa,maxiteropt,num,stopindex,r)
%% initialization algorithm
loss=[];
timerun=[0];
for i=1:num
    varze{i}=zeros(size(var{i},1),size(var{i},2));
end

cardinal=1e-6;
bts=[];
LK=zeros(1,num);
L=ones(1,num);
tk=1;
rate=[0];
varK=var;

wk=(tk-1)/(tk+1);


loss(1)=compute(var,num,ngmar);

returnloss=norm(ngmar,"fro");

t1=clock;



for i=1:maxiteropt

%% update parameters
fprintf("%d\n",i);
randsor=randperm(num);
randlen=randi([floor(num/2),num],1);
randsor=[randsor,randsample([1:num], randlen,true)];
varz=var;
LtempK=LK;
vartempK=varK;
varK=var;
LK=L;
[var,L,bt]=IBPGupdate(var,varze,vartempK,varz,num,ngmar,aa,randsor,L,LtempK,r,wk);
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
rate(i+1)=Res;
fprintf("cri:%d\n",Res);

stop=stopcheck(timerun,stopindex);
if(stop==1)
    fprintf("Number of terminations：%d\n",i);
    siter=i;
    pause(4);
    break;
end



tk=(1+sqrt(1+4*tk^2))/2;

wk=(tk-1)/(tk+1);

end




end













