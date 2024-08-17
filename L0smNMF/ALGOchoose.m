%% Input.
% var         : Initial matrix
% ngmar       : Decomposed matrix
% aa          : Maximum number of non-zero elements for each decomposition matrix
% The remaining parameters are explained the same as the 'main_Run_me' function

%% Output.
% vars:       : Decomposition matrix resulting from the final iterative result
% loss:       : Array of loss functions generated during iteration
% tr:         : Runtime array during iteration
% btss and atz: An array of extrapolated parameters produced by each algorithm during iteration,
%             : where btss represents the array of additional extrapolated parameters for algorithms, 
%             : and 'atz' represents the array of additional extrapolated parameters for IBPG and iPALM
%% An Accelerated Block Proximal method with Adaptive Momentum for Nonconvex and Nonsmooth Optimization

function [data,varss]=ALGOchoose(var,ngmar,aa,maxiteropt,bt,objs,num,flag,stopindex,r,t,btmax)
if(flag==0)
[vars,loss,tr,btss]=PALM(var,ngmar,aa,maxiteropt,num,stopindex,r);
varss{1}=vars;
lossdata{1}=loss;
trdata{1}=tr;
btz{1}=btss;

elseif(flag==1) 
[vars,loss,tr,btss,atz]=IPALM(var,ngmar,aa,maxiteropt,num,stopindex);
varss{1}=vars;
lossdata{1}=loss;
trdata{1}=tr;
btz{1}=btss;
data{4}=atz;

elseif(flag==2)
[vars,loss,tr,btss]=BPL(var,ngmar,aa,maxiteropt,num,stopindex,r);
varss{1}=vars;
lossdata{1}=loss;
trdata{1}=tr;
btz{1}=btss;


elseif(flag==3) 
[vars,loss,tr,btss,ats]=IBPG(var,ngmar,aa,maxiteropt,num,stopindex,r);
varss{1}=vars;
lossdata{1}=loss;
trdata{1}=tr;
btz{1}=btss;
atz{1}=ats;
data{4}=atz;



elseif(flag==4) 
[vars,loss,tr,btss]=TITAN(var,ngmar,aa,maxiteropt,num,stopindex,r);
varss{1}=vars;
lossdata{1}=loss;
trdata{1}=tr;
btz{1}=btss;


elseif(flag==5)
[vars,loss,tr,btss,ats]=TiBPALM(var,ngmar,aa,maxiteropt,num,stopindex,r); 
varss{1}=vars;
lossdata{1}=loss;
trdata{1}=tr;
btz{1}=btss;
atz{1}=ats;
data{4}=atz;

elseif(flag==6)
for i=1:(length(objs))
[vars,loss,tr,btss]=ABPL(var,ngmar,aa,maxiteropt,bt,objs(i),num,stopindex,r,t,btmax); %%ABPL-cyclic and ABPL+-cyclic
varss{i}=vars;
lossdata{i}=loss;
trdata{i}=tr;
btz{i}=btss;
end

elseif(flag==7)
for i=1:(length(objs))
[vars,loss,tr,btss]=ABPLRD(var,ngmar,aa,maxiteropt,bt,objs(i),num,stopindex,r,t,btmax); %%ABPL-random and ABPL+-random
varss{i}=vars;
lossdata{i}=loss;
trdata{i}=tr;
btz{i}=btss;
end


end








data{1}=lossdata;
data{2}=trdata;
data{3}=btz;

end



