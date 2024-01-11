%%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions
%%  This algorithm comes from the paper: "Inertial Proximal Alternating Linearized Minimization (iPALM) for Nonconvex and Nonsmooth Problems"
function [var,loss,timerun,bts,ats]=IPALM(var,ngmar,aa,maxiteropt,num,stopindex)
%% initialization algorithm
loss=[];
timerun=[0];
for i=1:num
    varze{i}=zeros(size(var{i},1),size(var{i},2));
end
bts=[];
varK=var;
loss(1)=compute(var,num,ngmar);
returnloss=norm(ngmar,"fro");


t1=clock;
for i=1:maxiteropt
fprintf("%d\n",i); 
%% This extrapolation parameter update strategy is the best choice in iPALM according to the numerical results.
at=(i-1)/(2*(i+2));
bt=(i-1)/(i+2);
%% Update parameters
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
fprintf("Rel：%d\n",Res);
stop=stopcheck(timerun,stopindex);
if(stop==1)
    fprintf("Number of terminations：%d\n",i);
    pause(4);
    break;
end
end

end


%% The parameters update function of iPALM
function var=IPALMupdate(var,varze,varK,varz,num,ngmar,aa,at,bt)
for j=1:num
V=gradipalm(var,varz,varK,at,bt,ngmar,j,num);
[var{j},~]=PROX(varze{j},V,aa(j));
end
end




