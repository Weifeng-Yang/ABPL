
% An Accelerated Block Proximal Framework with Adaptive Momentum for Nonconvex and Nonsmooth Optimization

function stop=stopcheck(timerun,stopindex)
    stop=0;
    
    if(stopindex==1)
        if(timerun(end)>50)
        stop=1;
        end
    elseif(stopindex==2)
        if(timerun(end)>100)
        stop=1;
        end
    end
end