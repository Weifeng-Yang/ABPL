clear
clc

% Parameter.
%   dimension : The number and dimension of matrices to be decomposed
%   index      : The data set to be used, when index=1, use  mycielskian12.mat dataset,
%                when index=2, use lp_ship12l dataset
%   bt        : Initial extrapolated parameter values for APGnc+, ABPL and ABPL+
%   btmax     : Maximum extrapolated parameter values for APGnc+, ABPL and ABPL+
%   t         : Decay and increment rates for extrapolated parameters for APGnc+, ABPL and ABPL+
%   r         : step factor
%   maxiteropt: maximum iteration alloted to the algorithm
%   trigger   : Whether to enable the indicator array of each algorithm, where
%               0∈trigger, enable the PALM algorithm
%               1∈trigger, enable the APGnc++ algorithm
%               2∈trigger, enable the iPALM algorithm
%               3∈trigger, enable the BPL algorithm
%               4∈trigger, enable the IBPG algorithm
%               5∈trigger, enable the ABPL+ algorithm
%               6∈trigger, enable the ABPL+-random algorithm 
%   objs      : Whether to enable ABPL+'s adaptive momentum
%               1∈objs, enable the adaptive momentum
%               2∈objs, enable the non-adaptive momentum
%   percent   ：The proportion of non-zero elements of each decomposition matrix
%   stopindex : The indicator of the stop condition, 
%               set the specific termination condition, see the 'stopcheck' function for details, 
%               The default termination condition is: each algorithm runs for 300 seconds
% Display
%   nonzero   ：The number of non-zero elements in each component.
%   cri       : The difference in the objective function value between two iterations.

%% Parameter settings
rng('shuffle')
warning('off');
dimension=[300,300];
index=2;
maxiteropt=60000;
r=1.5;
bt=0.6;
t=1.1;
btmax=0.9999;
trigger=[0,1,2,3,4,5,6];
objs=[1,2];
percent=0.3;
stopindex=2;






%% Algorithm iteration starts
[ngmar,dimension]=readfile(index,dimension);
for i=1:length(dimension)-1
    a(i)=dimension(i)*dimension(i+1)*percent;
end
aa=a;
num=length(dimension)-1;


for i=1:num
    var{i}=sprand(dimension(i),dimension(i+1),aa(i)/(2*dimension(i)*dimension(i+1)));
    var{i}=full(var{i});
end
for i=1:length(trigger)       
[datas{i},vars{i}]=ALGOchoose(var,ngmar,aa,maxiteropt,bt,objs,num,trigger(i),stopindex,r,t,btmax);
end


%% Drawing
plt=plotplt(datas,objs,trigger);







