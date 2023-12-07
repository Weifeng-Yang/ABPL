%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions

function [var,loss,timerun,bts,rate]=BPL(var,ngmar,aa,maxiteropt,num,N,stopindex,r)
%% initialization algorithm
loss=[];
timerun=[0];
for i=1:num
    varze{i}=zeros(size(ngmar,i),N);
end

LK=zeros(1,num);
L=ones(1,num);
tk=1;
cardinal=1e-8;
bts=[];
rate=[0];
wk=zeros(1,num);
varK=var;
for j=1:num
wk(j)=(tk-1)/(tk);
end

for j=1:num
    fprintf("nonzero:%d\n",sum(sum(var{j}~=0)));
end




returnloss=norm(ngmar);
loss(1)=compute(var,ngmar);


t1=clock;


for i=1:maxiteropt
%% update parameters
fprintf("%d\n",i);
if(i>1)
for j=1:num   
    wk(j)=min(wk(j),0.9999*(r-1)/(2*r+2)*sqrt(LK(j)/L(j)));
    vv{j}=var{j}+wk(j)*(var{j}-varK{j});
end
else
    vv=var;
end

varK=var;
LK=L;
var=vv;

randsor=randperm(num);
for k=1:num
    j=randsor(k);
    [V,L(j)]=gradbpl(var,ngmar,j,num,r);
    [var{j},~]=PROX(varze{j},var{j},V,aa(j));
end
loss(i+1)=compute(var,ngmar);


%% Judging whether to extrapolate
if(loss(i+1)>loss(i))
    var=varK;
    for k=1:num
    j=randsor(k);
    [V,L(j)]=gradbpl(var,ngmar,j,num,r);
    [var{j},~]=PROX(varze{j},var{j},V,aa(j));
    end
    loss(i+1)=compute(var,ngmar);
end



%% Check if termination condition is met
fprintf("BPL\n");
for j=1:num
    fprintf("nonzero:%d\n",nnz(var{j}~=0));
end
bts{i}=wk;
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

tk=(1+sqrt(1+4*tk^2))/2;
for j=1:num
wk(j)=(tk-1)/(tk);
end


end




end













