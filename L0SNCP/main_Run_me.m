clearvars
clc

% Parameter.
%   N         : Rank of Tensor CP Decomposition
%   index     : The data set to be used, when index=1, use  urban_traffic dataset,
%               when index=2, use PosteriorEmssion dataset
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
warning('off');
rng('shuffle')
N=50;
index=1;
bt=0.2;
btmax=0.9999;
t=1.3;
r=1.5;
maxiteropt=600000;
trigger=[3,4];
objs=[1,2];
percent=0.3;
stopindex=1;





%% Algorithm iteration starts
[ngmar,num]=readfile(index);

ngmar=tensor(ngmar);
for i=1:length(size(ngmar))
    a(i)=N*size(ngmar,i)*percent;
end
aa=a;
for i=1:num
    den=size(ngmar,i)*N;
    var{i}=sprand(size(ngmar,i),N,aa(i)/(2*den));
    var{i}=full(var{i});
end
for i=1:length(trigger)
[datas{i},vars{i}]=ALGOchoose(var,ngmar,aa,maxiteropt,bt,objs,num,N,trigger(i),stopindex,r,t,btmax);
end
datas{length(trigger)+1}=var;


%% Drawing
plt0=plotplt(datas,objs,trigger);








