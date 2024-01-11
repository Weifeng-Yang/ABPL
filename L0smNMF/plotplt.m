function plt0=plotplt(datas,objs,trigger)
color=["-o","-^","->","-d","-x","-p","-+"];
% color=["-o","-*","-+","-x","-d","-s","-p","-h"];
mess={};
plt0=0;
for i=1:length(trigger)
data=datas{i};
lossdata=data{1};
trdata=data{2};

if(~isempty(lossdata))
    mess=pltplt(lossdata,trdata,objs,trigger(i),mess,color(i));
end
end
end


function mess=pltplt(lossdata,trdata,objs,trigger,mess,color)
figure(1)
ss=10;
for i=1:length(lossdata)
    maker_idx = 1:ss:length(lossdata{i});
    semilogy(trdata{i},lossdata{i},color,'linewidth',2,'MarkerIndices',maker_idx);
    xlabel('Time (seconds)','FontSize',25);
    ylabel('Objective funciton value','FontSize',25);
    if(trigger==0)
     mes{i}='PALM';
    elseif(trigger==1)
    
    mes{i}='iPALM';
    
    elseif(trigger==2)
    mes{i}='BPL';
    elseif(trigger==3)
    mes{i}='IBPG';
   elseif(trigger==4)
    mes{i}='TITAN';
    elseif(trigger==5)
        if(objs(i)==1)
            mes{i}='ABPL$^{+}$-cyclic';
        elseif(objs(i)==2)
            mes{i}='ABPL-cyclic';
        end
    elseif(trigger==6)
       if(objs(i)==1)
            mes{i}='ABPL$^{+}$-random';
        elseif(objs(i)==2)
            mes{i}='ABPL-random';
       end   
    end
    hold on;
end
ls=length(mess);
for j=1:length(mes)
    mess{j+ls}=mes{j};
end
h=legend(mess,'Interpreter','latex');
set(gca,'FontSize',50);
set(h,'FontSize',28)
xlabel('Time (seconds)','FontSize',50);
ylabel('Objective funciton value','FontSize',50);
prettyAxes().gbase2()
end




