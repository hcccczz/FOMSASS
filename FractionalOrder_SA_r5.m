function [Indexes,Indexes0,Indexes1,Indexes2,Indexes3,Indexes4,Popul0,Popul1,Popul2,Popul3,Popul4,Archive0,Archive1,Archive2,Archive3,Archive4,Popul_fo,Popul_best,Archive_fo]=FractionalOrder_SA_r5(Indexes,Indexes0,Indexes1,Indexes2,Indexes3,Indexes4,Popul0,Popul1,Popul2,Popul3,Popul4,Archive0,Archive1,Archive2,Archive3,Archive4,pop,Archive,Popul_fo,Popul_best,Archive_fo,fo_rate,generation)
%%更新分数阶信息
%% 这是保留五层记忆
ArchiveSize=size(Archive,1);
[popsize,D]=size(pop);
% if(size(Archive5,1)<ArchiveSize)
%     temp=zeros(ArchiveSize-size(Archive5,1),D);
%     Archive5=[Archive5;temp];
% 
% end
if(size(Archive4,1)<ArchiveSize)
    temp=zeros(ArchiveSize-size(Archive4,1),D);
    Archive4=[Archive4;temp];

end
if(size(Archive3,1)<ArchiveSize)
    temp=zeros(ArchiveSize-size(Archive3,1),D);
    Archive3=[Archive3;temp];

end
if(size(Archive2,1)<ArchiveSize)
    temp=zeros(ArchiveSize-size(Archive2,1),D);
    Archive2=[Archive2;temp];

end
if(size(Archive1,1)<ArchiveSize)
    temp=zeros(ArchiveSize-size(Archive1,1),D);
    Archive1=[Archive1;temp];

end
if(generation==1)
    Popul_fo=pop;
    Popul_best=pop(Indexes,:);
    Popul1=pop;
    Archive_fo=Archive;
    Archive1=Archive;
    Indexes1=Indexes;
elseif(generation==2)
    Popul_fo=(1/gamma(2)*fo_rate)*pop+1/gamma(3)*fo_rate*(1-fo_rate)*Popul1(1:popsize,:);
    Popul_best=(1/gamma(2)*fo_rate)*pop(Indexes,:)+1/gamma(3)*fo_rate*(1-fo_rate)*Popul1(Indexes1(1:popsize),:);
    Popul2=Popul1;
    Popul1=pop;
    Archive_fo=(1/gamma(2)*fo_rate)*Archive(1:ArchiveSize,:)+1/gamma(3)*fo_rate*(1-fo_rate)*Archive1(1:ArchiveSize,:);
    Archive2=Archive1;
    Archive1=Archive;
    Indexes2=Indexes1;
    Indexes1=Indexes;
elseif(generation==3)
    Popul_fo=(1/gamma(2)*fo_rate)*pop+1/gamma(3)*fo_rate*(1-fo_rate)*Popul1(1:popsize,:)+1/gamma(4)*fo_rate*(1-fo_rate)*(2-fo_rate)*Popul2(1:popsize,:);
    Popul_best=(1/gamma(2)*fo_rate)*pop(Indexes,:)+1/gamma(3)*fo_rate*(1-fo_rate)*Popul1(Indexes1(1:popsize),:)+1/gamma(4)*fo_rate*(1-fo_rate)*(2-fo_rate)*Popul2(Indexes2(1:popsize),:);
    Popul3=Popul2;
    Popul2=Popul1;
    Popul1=pop;
    Archive_fo=(1/gamma(2)*fo_rate)*Archive(1:ArchiveSize,:)+1/gamma(3)*fo_rate*(1-fo_rate)*Archive1(1:ArchiveSize,:)+1/gamma(4)*fo_rate*(1-fo_rate)*(2-fo_rate)*Archive2(1:ArchiveSize,:);
    Archive3=Archive2;
    Archive2=Archive1;
    Archive1=Archive;
    Indexes3=Indexes2;
    Indexes2=Indexes1;
    Indexes1=Indexes;
elseif(generation==4)
    Popul_fo=(1/gamma(2)*fo_rate)*pop+1/gamma(3)*fo_rate*(1-fo_rate)*Popul1(1:popsize,:)+1/gamma(4)*fo_rate*(1-fo_rate)*(2-fo_rate)*Popul2(1:popsize,:)+...
        1/gamma(5)*fo_rate*(1-fo_rate)*(2-fo_rate)*(3-fo_rate)*Popul3(1:popsize,:);
    Popul_best=(1/gamma(2)*fo_rate)*pop(Indexes,:)+1/gamma(3)*fo_rate*(1-fo_rate)*Popul1(Indexes1(1:popsize),:)+1/gamma(4)*fo_rate*(1-fo_rate)*(2-fo_rate)*Popul2(Indexes2(1:popsize),:)+...
        1/gamma(5)*fo_rate*(1-fo_rate)*(2-fo_rate)*(3-fo_rate)*Popul3(Indexes3(1:popsize),:);
    Popul4=Popul3;
    Popul3=Popul2;
    Popul2=Popul1;
    Popul1=pop;
    Archive_fo=(1/gamma(2)*fo_rate)*Archive+1/gamma(3)*fo_rate*(1-fo_rate)*Archive1(1:ArchiveSize,:)+1/gamma(4)*fo_rate*(1-fo_rate)*(2-fo_rate)*Archive2(1:ArchiveSize,:)+...
        1/gamma(5)*fo_rate*(1-fo_rate)*(2-fo_rate)*(3-fo_rate)*Archive3(1:ArchiveSize,:);
    Archive4=Archive3;
    Archive3=Archive2;
    Archive2=Archive1;
    Archive1=Archive;
    Indexes4=Indexes3;
    Indexes3=Indexes2;
    Indexes2=Indexes1;
    Indexes1=Indexes;
else
    Popul_fo=(1/gamma(2)*fo_rate)*pop+1/gamma(3)*fo_rate*(1-fo_rate)*Popul1(1:popsize,:)+1/gamma(4)*fo_rate*(1-fo_rate)*(2-fo_rate)*Popul2(1:popsize,:)+...
        1/gamma(5)*fo_rate*(1-fo_rate)*(2-fo_rate)*(3-fo_rate)*Popul3(1:popsize,:)+...
        1/gamma(6)*fo_rate*(1-fo_rate)*(2-fo_rate)*(3-fo_rate)*(4-fo_rate)*Popul4(1:popsize,:);
    Popul_best=(1/gamma(2)*fo_rate)*pop(Indexes,:)+1/gamma(3)*fo_rate*(1-fo_rate)*Popul1(Indexes1(1:popsize),:)+1/gamma(4)*fo_rate*(1-fo_rate)*(2-fo_rate)*Popul2(Indexes2(1:popsize),:)+...
        1/gamma(5)*fo_rate*(1-fo_rate)*(2-fo_rate)*(3-fo_rate)*Popul3(Indexes3(1:popsize),:)+...
        1/gamma(6)*fo_rate*(1-fo_rate)*(2-fo_rate)*(3-fo_rate)*(4-fo_rate)*Popul4(Indexes4(1:popsize),:);
    % Popul5=Popul4;
    Popul4=Popul3;
    Popul3=Popul2;
    Popul2=Popul1;
    Popul1=pop;
    Archive_fo=(1/gamma(2)*fo_rate)*Archive+1/gamma(3)*fo_rate*(1-fo_rate)*Archive1(1:ArchiveSize,:)+1/gamma(4)*fo_rate*(1-fo_rate)*(2-fo_rate)*Archive2(1:ArchiveSize,:)+...
        1/gamma(5)*fo_rate*(1-fo_rate)*(2-fo_rate)*(3-fo_rate)*Archive3(1:ArchiveSize,:)+...
        1/gamma(6)*fo_rate*(1-fo_rate)*(2-fo_rate)*(3-fo_rate)*(4-fo_rate)*Archive4(1:ArchiveSize,:);
    % Archive5=Archive4;
    Archive4=Archive3;
    Archive3=Archive2;
    Archive2=Archive1;
    Archive1=Archive;
    % Indexes5=Indexes4;
    Indexes4=Indexes3;
    Indexes3=Indexes2;
    Indexes2=Indexes1;
    Indexes1=Indexes;
% else
%     Popul_fo=(1/gamma(2)*fo_rate)*pop+1/gamma(3)*fo_rate*(1-fo_rate)*Popul1(1:popsize,:)+1/gamma(4)*fo_rate*(1-fo_rate)*(2-fo_rate)*Popul2(1:popsize,:)+...
%         1/gamma(5)*fo_rate*(1-fo_rate)*(2-fo_rate)*(3-fo_rate)*Popul3(1:popsize,:)+...
%         1/gamma(6)*fo_rate*(1-fo_rate)*(2-fo_rate)*(3-fo_rate)*(4-fo_rate)*Popul4(1:popsize,:)+...
%         1/gamma(7)*fo_rate*(1-fo_rate)*(2-fo_rate)*(3-fo_rate)*(4-fo_rate)*(5-fo_rate)*Popul5(1:popsize,:);
%     Popul_best=(1/gamma(2)*fo_rate)*pop(Indexes,:)+1/gamma(3)*fo_rate*(1-fo_rate)*Popul1(Indexes1(1:popsize),:)+1/gamma(4)*fo_rate*(1-fo_rate)*(2-fo_rate)*Popul2(Indexes2(1:popsize),:)+...
%         1/gamma(5)*fo_rate*(1-fo_rate)*(2-fo_rate)*(3-fo_rate)*Popul3(Indexes3(1:popsize),:)+...
%         1/gamma(6)*fo_rate*(1-fo_rate)*(2-fo_rate)*(3-fo_rate)*(4-fo_rate)*Popul4(Indexes4(1:popsize),:)+...
%         1/gamma(7)*fo_rate*(1-fo_rate)*(2-fo_rate)*(3-fo_rate)*(4-fo_rate)*(5-fo_rate)*Popul5(Indexes5(1:popsize),:);
%     Popul5=Popul4;
%     Popul4=Popul3;
%     Popul3=Popul2;
%     Popul2=Popul1;
%     Popul1=pop;
%     Archive_fo=(1/gamma(2)*fo_rate)*Archive+1/gamma(3)*fo_rate*(1-fo_rate)*Archive1(1:ArchiveSize,:)+1/gamma(4)*fo_rate*(1-fo_rate)*(2-fo_rate)*Archive2(1:ArchiveSize,:)+...
%         1/gamma(5)*fo_rate*(1-fo_rate)*(2-fo_rate)*(3-fo_rate)*Archive3(1:ArchiveSize,:)+...
%         1/gamma(6)*fo_rate*(1-fo_rate)*(2-fo_rate)*(3-fo_rate)*(4-fo_rate)*Archive4(1:ArchiveSize,:)+...
%         1/gamma(7)*fo_rate*(1-fo_rate)*(2-fo_rate)*(3-fo_rate)*(4-fo_rate)*(5-fo_rate)*Archive5(1:ArchiveSize,:);
%     Archive5=Archive4;
%     Archive4=Archive3;
%     Archive3=Archive2;
%     Archive2=Archive1;
%     Archive1=Archive;
%     Indexes5=Indexes4;
%     Indexes4=Indexes3;
%     Indexes3=Indexes2;
%     Indexes2=Indexes1;
%     Indexes1=Indexes;

    





end


end

