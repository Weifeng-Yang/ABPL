function check=plotplt(datas,objs,trigger)
color=["-o","-*","-+","-x","-p","-d","->"];
mess={};
check=0;
num=length(datas{end});
for i=1:length(trigger)
data=datas{i};
lossdata=data{1};
trdata=data{2};
rates=data{3};
bts=data{4};

if(~isempty(lossdata))
    mess=pltplt(lossdata,trdata,objs,trigger(i),mess,color(i));
end
end
end

function mess=pltplt(lossdata,trdata,objs,trigger,mess,color)
figure(1)
% color=["-o","-*","-+","-x","-d","-s","-p","-h"];

ss=20;
for i=1:length(lossdata)
    maker_idx = 1:ss:length(lossdata{i});
    semilogy(trdata{i},lossdata{i},color,'linewidth',0.9,'MarkerIndices',maker_idx);
    xlabel('iter(count)','FontSize',30);
    ylabel('Objective funciton value','FontSize',30);
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
set(h,'FontSize',12)
% prettyAxes().dark()

figure(2)
for i=1:length(lossdata)
        maker_idx = 1:ss:length(lossdata{i});
    semilogy(trdata{i},lossdata{i},color,'linewidth',0.9,'MarkerIndices',maker_idx);
    xlabel('Time(S)','FontSize',25);
    ylabel('Objective funciton value','FontSize',25);
    hold on;
end
h=legend(mess,'Interpreter','latex');
set(h,'FontSize',11)
% figure(3)
% for i=1:length(lossdata)
%     plot(trdata{i},rates{i},color);
%     xlabel(['rumtime(S),dim=',num2str(dim)]);
%     ylabel('RelErr');
%     hold on;
% end
h=legend(mess,'Interpreter','latex');
set(h,'FontSize',11)
    plt0=0;
end

