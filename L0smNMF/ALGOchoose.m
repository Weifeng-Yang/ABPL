
% Input.
% var         : initial matrix
% ngmar       : decomposed matrix
% aa          : Maximum number of non-zero elements for each decomposition matrix
% The remaining parameters are explained the same as the 'main_Run_me' function

% Output.
% vars:       : Decomposition matrix resulting from the final iterative result
% loss:       : Array of loss functions generated during iteration
% tr:         :Runtime array during iteration
% btss and atz:An array of extrapolated parameters produced by each algorithm during iteration
%             :where 'atz' represents the array of additional extrapolated parameters for IBPG and iPALM
% rate        :Decomposed residuals; see Eq. (24) of the paper 
% An Accelerated Block Proximal Framework with Adaptive Momentum for Nonconvex and Nonsmooth Optimization

function [data,varss]=ALGOchoose(var,ngmar,aa,maxiteropt,bt,objs,num,flag,stopindex,r,t,btmax)
if(flag==0)
[vars,loss,tr,btss,rate]=PALM(var,ngmar,aa,maxiteropt,num,stopindex,r);
varss{1}=vars;
lossdata{1}=loss;
trdata{1}=tr;
btz{1}=btss;
rates{1}=rate;

elseif(flag==1) 
[vars,loss,tr,btss,atz,rate]=IPALM(var,ngmar,aa,maxiteropt,bt,num,stopindex);
varss{1}=vars;
lossdata{1}=loss;
trdata{1}=tr;
btz{1}=btss;
rates{1}=rate;
data{5}=atz;

elseif(flag==2)
[vars,loss,tr,btss,rate]=BPL(var,ngmar,aa,maxiteropt,num,stopindex,r);
varss{1}=vars;
lossdata{1}=loss;
trdata{1}=tr;
btz{1}=btss;
rates{1}=rate;


elseif(flag==3) 
[vars,loss,tr,btss,ats,rate]=IBPG(var,ngmar,aa,maxiteropt,num,stopindex,r);
varss{1}=vars;
lossdata{1}=loss;
trdata{1}=tr;
btz{1}=btss;
rates{1}=rate;
atz{1}=ats;
data{5}=atz;



elseif(flag==4) 
[vars,loss,tr,btss,rate]=TITAN(var,ngmar,aa,maxiteropt,num,stopindex,r);
varss{1}=vars;
lossdata{1}=loss;
trdata{1}=tr;
btz{1}=btss;
rates{1}=rate;


elseif(flag==5)
for i=1:(length(objs))
[vars,loss,tr,btss,rate]=ABPL(var,ngmar,aa,maxiteropt,bt,objs(i),num,stopindex,r,t,btmax); %%ABPL-cyclic and ABPL+-cyclic
varss{i}=vars;
lossdata{i}=loss;
trdata{i}=tr;
btz{i}=btss;
rates{i}=rate;
data{5}=objs;
end

elseif(flag==6)
for i=1:(length(objs))
[vars,loss,tr,btss,rate]=ABPLRD(var,ngmar,aa,maxiteropt,bt,objs(i),num,stopindex,r,t,btmax); %%ABPL-random and ABPL+-random
varss{i}=vars;
lossdata{i}=loss;
trdata{i}=tr;
btz{i}=btss;
rates{i}=rate;
data{5}=objs;
end


end





data{1}=lossdata;
data{2}=trdata;
data{3}=rates;
data{4}=btz;





end