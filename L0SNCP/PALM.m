%%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions
%%  This algorithm comes from the paper: "Proximal alternating linearized minimization for nonconvex and nonsmooth problems"
function [var,loss,timerun,bts]=PALM(var,ngmar,aa,maxiteropt,num,N,stopindex,r)
%% initialization algorithm
loss=[];
timerun=[0];
for i=1:num
    varze{i}=zeros(size(ngmar,i),N);
end
bts=0;
for j=1:num
    fprintf("nonzero:%d\n",sum(sum(var{j}~=0)));
end
returnloss=norm(ngmar);
loss(1)=compute(var,ngmar);
t1=clock;


for i=1:maxiteropt
%% update parameters
fprintf("%d\n",i);
varK=var;
for j=1:num
V=grad(var,ngmar,j,num,r);
var{j}=PROX(varze{j},V,aa(j));
end
loss(i+1)=compute(var,ngmar);


%% Check if termination condition is met
fprintf("PALM\n");
for j=1:num
    fprintf("nonzero:%d\n",nnz(var{j}~=0));
end
t2=clock;
timerun(i+1)=etime(t2,t1);
ReeK=loss(i)/returnloss;
Ree=loss(i+1)/returnloss;
Res=abs(Ree-ReeK);
fprintf("Rel：%d\n",Res);
stop=stopcheck(timerun,stopindex);
if(stop==1)
    fprintf("Number of terminations：%d\n",i);
    pause(4);
    break;
end
end

end















