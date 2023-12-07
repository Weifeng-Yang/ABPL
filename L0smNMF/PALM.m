%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions

function [var,loss,timerun,bts,rate]=PALM(var,ngmar,aa,maxiteropt,num,stopindex,r)

%% initialization algorithm
loss=[];
timerun=[0];
rate=[0];
for i=1:num
    varze{i}=zeros(size(var{i},1),size(var{i},2));
end

bts=0;
varK=var;


loss(1)=compute(var,num,ngmar);

cardinal=1e-6;
returnloss=norm(ngmar,"fro");


t1=clock;
for i=1:maxiteropt
fprintf("%d\n",i);



%% update parameters
var=paraupdate(var,varze,num,ngmar,aa,r);
loss(i+1)=compute(var,num,ngmar);



t2=clock;
timerun(i+1)=etime(t2,t1);





%% Check if termination condition is met
fprintf("PALM\n");
for j=1:num
    fprintf("nonzero:%d\n",nnz(var{j}));
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




