clc;
clear all;
pp = [10,30,50,100];
%% Loop for different dimension problems
for jj = 1:4
format long;
format compact;
FF = zeros(51,30);
D = pp(jj); %问题维度
Max_NFES = 10000 * D;
val_2_reach = 10^(-8);
Max_Region = 100.0;
Min_Region = -100.0;
LU = [-100 * ones(1, D); 100 * ones(1, D)]; 
fhd=@cec17_func;
pb = 0.4;
ps = 0.5;
num_prbs = 30;
runs = 51;
run_funcvals = [];
RecordFEc_iactor = ...
    [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, ...
    0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
progress = numel(RecordFEc_iactor);
allerrorvals = zeros(progress, runs, num_prbs);
result=zeros(num_prbs,5);
fprintf('Running FOSMSASS on D= %d\n', D)
%% Loop for different problems
for func = 30 : num_prbs
    Optm = func * 100.0;
    outcome = [];
    fprintf('\n-------------------------------------------------------\n')
    fprintf('Function = %d, Dimension size = %d\n', func, D)
    

    for run_id = 1 : runs
        rand('seed',run_id*sum(100*clock));
        run_funcvals = [];
        col=1;
        generation=0;

        PbestRate = 0.11;
        popsize = 18*D;
        % popsize = 100;
        % rd = 0.8;
        rd = 0.5;

        % c=0.3;
        c=0.7;

        Ar=1.3;

         Ms = 6;  %档案数是6这里
         nA=0;
         
        % popsize=floor( 25*log(D)*sqrt(D));
        SEL = round(ps*popsize);
        Max_Pop_Size = popsize;
        Min_Pop_Size=4.0;
        G_Max = ceil(log(Min_Pop_Size/Max_Pop_Size)/log(1+(Min_Pop_Size-Max_Pop_Size)/Max_NFES));  %最大代数
        NFES = 0;
        ProblemSize=D;
        %% 开始
        % fo_rate=Parameter.fo_rate;
        % fo% fo_rate=0.99;
        % fo_rate=0.999;

        % fo_rate=0.9;
        % fo_rate=0.7;
        % fo_rate=0.5;
        % fo_rate=0.1;

        
        
        convergence=[];
        %% 这是分数阶的参数
        Indexes = [];
        Indexes0 =  []; %分数阶微积分相关索引
        Indexes1 = [];
        Indexes2 = [];
        Indexes3 =  [];
        Indexes4 = [];
        Indexes5 = [];
        Popul0=[];
        Popul1=[];
        Popul2=[];
        Popul3=[];
        Popul4=[];
        Popul5=[];
        Popul_fo=[]; %分数阶导数所用
        Popul_best=[];%分数阶导数所用
        Archive0=[];
        Archive1=[];
        Archive2=[];
        Archive3=[];
        Archive4=[];
        Archive5=[];
        Archive_fo=[]; %分数阶导数所用


        % %% SSAP 自适应的两个种群的最优值
        % i1BestFit = 1e+30; %原文中的F(xr)，前一部分群体中最好的适应度
        % i2BestFit = 1e+30;    %原文中的F(xb),后一部分群体中最好的适应度
        % l = 2;
        % lList=[];
        % wn1 = 0;
        % wn2 = 0;
        % wn3 = 0;
        convergence=[];
        % % 自适应的参数
        % Msize=6; %选取1 6 50 试一下  试一下6  10  20 50 100
        % M_p=zeros(Msize,3);  %这是SSM中的M_p矩阵，第一列存放的是上一轮计算得到的概率因子rk+1  第二列存放的是优化的适应度均值
        % M_p=0.5*ones(Msize,1);
        % M_p(:,1)=0.5; % 初始情况设置人口分配比设置为0.5
        % M_iterator=1;

       
        
       %% Initialize the main population
        PopOld = repmat(LU(1, :), popsize, 1) + rand(popsize, ProblemSize) .* (repmat(LU(2, :) - LU(1, :), popsize, 1));
        % PopOld2 = repmat(LU(1, :), popsize, 1) + rand(popsize, ProblemSize) .* (repmat(LU(2, :) - LU(1, :), popsize, 1));
        Pop = PopOld; % the old Population becomes the current Population
        fitness = feval(fhd,Pop',func); %计算适应度，返回的是一个行向量，需要变成列向量
        fitness = fitness';
        % fitness_temp=fitness;
        start_record=min(fitness);
        bc_i_fit_var = 1e+30; 
        bc_i_index = 0; 
        bc_i_solution = zeros(1, D);
        for i = 1 : popsize
            NFES = NFES + 1;
            if fitness(i) < bc_i_fit_var
                bc_i_fit_var = fitness(i); %最好的适应度
                bc_i_solution = Pop(i, :);  %最好的个体的位置信息
                bc_i_index = i;  %最好的索引
            end
            if NFES > Max_NFES; break; end
        end
        
        MemCI = c .* ones(Ms, 1); %生成存放c记忆的列向量
        MemRD = rd .* ones(Ms, 1); %生成存放rank记忆的列向量
        MemI = 1;  %迭代器位置
        Arch.NP = Ar * popsize; 
        Arch.pop = zeros(0, ProblemSize); 
        Arch.OF = zeros(0, 1);
        %==========================================================================
        %==========================================================================
while NFES < Max_NFES
    % gg = gg+1;
    generation = generation+1;
    % PbestRate=0.085+0.085*NFES/Max_NFES;
    PbestRate=0.11;
    fo_rate=0.9+0.09*(NFES/Max_NFES); %0.9-0.99
    %% calculation of z for every individuals
    n = floor(popsize/2);  %选取一半个体，前一半个体采用随机的变异方式，后一半采用Pbest的变异方式
    % PopOld2 = PopOld2(randperm(popsize),:);
    
    %% 计算种群多样性
    % diver_save(generation,run_id)=Diversity2(PopOld);



    % l=mean(M_p);
    % L_save(generation,run_id)=l;
    % lList(generation)=l;
    % n=ceil(popsize*l);

    Pop = PopOld; 
    [~, Indexes] = sort(fitness, 'ascend');
    Indexes=Indexes.';
    ks = Indexes(1:n); %前一半个体的索引 better
    worst_index=Indexes(n+1:end); % 后一半个体的索引 worse
    MemCI(Ms)=0.9;
    MemRD(Ms)=0.9;

    mem_rand_index = ceil(Ms * rand(popsize, 1)); %获得每个个体的记忆的随机索引
    MUci = MemCI(mem_rand_index); %获得每个个体选取的miuCI值
    MUrd = MemRD(mem_rand_index); %获得每个个体的miurand值
    rd = normrnd(MUrd, 0.1);      %使用正态分布获得了rand值，rand为一个列向量，为popsizex1
    Term_Pos = find(MUrd == -1);  %检查有没有小于0的情况
    rd(Term_Pos) = 0;             %小于0的情况设置为0
    rd = min(rd, 1);              %将rd限制在0-1内
    rd = max(rd, 0);
    ci = MUci + 0.1 * tan(pi * (rand(popsize, 1) - 0.5));  %ci步长采用了柯西分布
    Pos = find(ci <= 0); %如果ci小于就重新进行柯西分布计算ci
    while ~ isempty(Pos)
	     ci(Pos) = MUci(Pos) + 0.1 * tan(pi * (rand(length(Pos), 1) - 0.5));
	     Pos = find(ci <= 0);
    end
    ci = min(ci, 1);  %防止ci越界，隔断在1
    %% 对参数的限制 ，ci rd
    if NFES<0.6*Max_NFES
        ci(ci>0.7) = 0.7;
    end
    
    a=1;
    if NFES<0.25*Max_NFES
        rd(rd<0.7) = 0.7;
    elseif NFES<0.5*Max_NFES
        rd(rd<0.6) = 0.6;
    end

    r0 = 1 : popsize;

    %% 分数阶模块
    Archive=Arch.pop;
    %%下面生成随机数
    [Indexes,Indexes0,Indexes1,Indexes2,Indexes3,Indexes4,Indexes5,Popul0,Popul1,Popul2,Popul3,Popul4,Popul5,Archive0,Archive1,Archive2,Archive3,Archive4,Archive5,Popul_fo,Popul_best,Archive_fo]=...
        FractionalOrder_SA(Indexes,Indexes0,Indexes1,Indexes2,Indexes3,Indexes4,Indexes5,Popul0,Popul1,Popul2,Popul3,Popul4,Popul5,Archive0,Archive1,Archive2,Archive3,Archive4,Archive5,Pop,Archive,Popul_fo,Popul_best,Archive_fo,fo_rate,generation);




    % PopAll = [Pop; Arch.pop];  %档案个体
    % PopAll = [Popul_fo; Archive_fo];
    % [r1, r2, r3] = genR1R2R3(popsize, size(PopAll, 1), r0); %r2是从档案和种群的并集中选出的索引，r1，r2为随机索引
    % pNP = max(round(PbestRate * popsize), 2);  %pbest的个数，并且最小为2
    % randindex = ceil(rand(1, popsize) .* pNP);  %获得pbest的索引
    % randindex = max(1, randindex);                %防止出现0索引的情况
    % pbest = Pop(SortedIndex(randindex), :);     %获得pbest的个体
    % tep = pbest;tep(ks,:)=Pop(r3(ks),:);PopA = PopOld2(r1,:);PopA(ks,:)=PopAll(r2(ks),:); %前一半变异个体采用随机变异，tep即为p个体，前一半采用的p个体为随机，后一半个体采用的p个体为pbest，并且前一半个体的popA要使用档案

   
    % pNP = max(round(PbestRate * popsize), 2);  %pbest的个数，并且最小为2
    % randindex = ceil(rand(1, popsize) .* pNP);  %获得pbest的索引
    % pbest = Popul_best(randindex, :); %% randomly choose one of the top 100p% solutions 获得pbest个体的位置 ,Popul_best是已经排好序号
    pNP = max(round(PbestRate * popsize), 2); %% choose at least two best solutions
    randindex = ceil(rand(1, popsize) .* pNP); %% select from [1, 2, 3, ..., pNP] 从p个最好个精英个体选取pbest
    randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0 避免选取到0
    pos=Indexes(randindex)==r0; %这里这样用会出错，因为indexes为列向量，而r0为行向量
    while(sum(pos)~=0&&NFES/Max_NFES<0.5)
        randindex(pos)= ceil(rand(1, sum(pos)) .* pNP);
        randindex(pos) = max(1, randindex(pos));
        pos=Indexes(randindex)==r0;
    end
    pbest = Popul_best(randindex, :); %% randomly choose one of the top 100p% solutions 获得pbest个体的位置
    
     % [r1, r2, r3] = genR1R2R3(popsize, size(PopAll, 1), r0); 

    RSP_template=1:popsize;
    RSP=3.0.*(popsize-RSP_template); %rsp向量
    prob = RSP / sum(RSP); % 这是按顺序的排序，但是种群是没有进行排序的
    % prob(Indexes)=prob; %置放对应位置的选择顺序
    %% 这个版本的
    % [r1, r2,r3]=gnR1R2_prob2(r0,prob,Indexes);
    % [r1, r2, r3] = genR1R2R3_r1RSP(popsize, size(PopAll, 1), r0,prob); %r2是从档案和种群的并集中选出的索引，只对r1进行了压力选择
    [r1, r2,r3]=gnR1R2_prob2(r0,prob,Indexes); %请注意r0和Indexes的维度必须相同，即同时为行向量
    
    
    % archivepop=Archive_fo;
    % archiveSize=size(archivepop,1);

    archivepop=Archive_fo;
    archiveSize=size(archivepop,1);
    archiveFitness=Arch.OF;

    tep = pbest;
    tep(ks,:)=Popul_fo(r3(ks),:);
    zi=zeros(size(Pop));
    popR2=zeros(size(Pop)); %这是排完序后的r2
    fitnessR2=zeros(popsize,1); %需要保证适应度都为列向量
    useArchive=rand(1,popsize)<archiveSize/(archiveSize+popsize); 
    notuseArch=~useArchive;
    if(archiveSize>0)
        r2_arch=floor(rand(1,popsize).*archiveSize)+1;
        popR2(useArchive,:)=archivepop(r2_arch(useArchive),:);
        fitnessR2(useArchive,:)=archiveFitness(r2_arch(useArchive),:);
    end
    notuseArch=~useArchive;
    popR2(notuseArch,:)=Popul_fo(r2(notuseArch),:);
    fitnessR2(notuseArch,:)=fitness(r2(notuseArch),:);

    rand_temp=Popul_fo(r3(ks),:);
    rand_fit=fitness(r3(ks),:);
    r1_temp=Popul_fo(r1(ks),:);
    r1_fit=fitness(r1(ks),:);
    r2_temp=popR2(ks,:);
    r2_fit=fitnessR2(ks,:);
    fit_temp=[r1_fit r2_fit];
    [~,bestfit_index]=sort(fit_temp,2); 
    % rand_final=zeros(size(rand_temp));
    r1_final=zeros(size(r1_temp));
    r2_final=zeros(size(r2_temp));

    %% 先是r1temp的放置
    pos1=bestfit_index==1;
    % rand_final(pos1(:,1),:)=rand_temp(pos1(:,1),:);
    r1_final(pos1(:,1),:)=r1_temp(pos1(:,1),:);
    r2_final(pos1(:,2),:)=r1_temp(pos1(:,2),:);

    %%  然后是r2temp的放置，同理 ，2代表r1
    pos2=bestfit_index==2;
    % rand_final(pos2(:,1),:)=r1_temp(pos2(:,1),:);
    r1_final(pos2(:,1),:)=r2_temp(pos2(:,1),:);
    r2_final(pos2(:,2),:)=r2_temp(pos2(:,2),:);


    %%优势种群个体的三个zi分量排序完成
    zi(ks,:)=rand_temp-Popul_fo(r0(ks),:)+r1_final-r2_final;

    %% 然后是worst个体的zi
    r1_temp=Popul_fo(r1(worst_index),:);
    r1_fit=fitness(r1(worst_index),:);
    r2_temp=popR2(worst_index,:);
    r2_fit=fitnessR2(worst_index,:);
    fit_temp=[r1_fit r2_fit];
    [~,worstfit_index]=sort(fit_temp,2); %按行排序 ，worstfit_index表示即索引 ，第一列即为最好适应度的索引，第二列为最差适应度的索引；

    r1_final=zeros(size(r1_temp));
    r2_final=zeros(size(r2_temp));

    %% 先完成 r1_temp的放置，1代表r1
    pos1=worstfit_index==1;
    r1_final(pos1(:,1),:)=r1_temp(pos1(:,1),:);
    r2_final(pos1(:,2),:)=r1_temp(pos1(:,2),:);

    %% 再完成r2_temp的放置，2代表r2
    pos2=worstfit_index==2;
    r1_final(pos2(:,1),:)=r2_temp(pos2(:,1),:);
    r2_final(pos2(:,2),:)=r2_temp(pos2(:,2),:);

    %%worst个体的r1 r2分量排序完成
    zi(worst_index,:)=tep(worst_index,:)-Popul_fo(r0(worst_index),:)+r1_final-r2_final;

    %% zi生成完毕，下面开始投影操作

    % tep = pbest;
    % tep(ks,:)=Popul_fo(r3(ks),:);
    % % PopA = pop(order2,:);
    % % PopA(ks,:)=PopAll(r2(ks),:); 
    % 
    % useArchive=rand(1,popsize)<archiveSize/(archiveSize+popsize);
    % %%计算zi
    % if(archiveSize>0)
    %     r2_arch=floor(rand(1,popsize).*archiveSize)+1;
    %     zi(useArchive,:)=tep(useArchive,:) - Popul_fo(r0(useArchive),:) + Popul_fo(r1(useArchive), :) - archivepop(r2_arch(useArchive),:);
    % end
    % notuseArch=~useArchive;
    % % zi = tep - pop(r0,:) + pop(r1, :) - PopAll(r2,:); 
    % zi(notuseArch,:)=tep(notuseArch,:) - Popul_fo(r0(notuseArch),:) + Popul_fo(r1(notuseArch), :) - Popul_fo(r2(notuseArch),:);
    %% calculation of Orthogonal matrix  计算随机的正交矩阵
    A = RandOrthMat(ProblemSize,1e-12);  %获得正交矩阵A DXD
    %% calculation of yi = Pop + ci.A.diag(bi).A'zi in parallel
      zi = zi*A; 
      Ur = zeros(popsize, ProblemSize);  
      J = (mod(floor(rand(popsize, 1)*ProblemSize), ProblemSize))*popsize + (1:popsize)';
      bi = rand(popsize, ProblemSize) < rd(:, ones(1, ProblemSize));
      Ur(J) = zi(J);
      Ur(bi) = zi(bi);
      yi = Pop(r0,:) + ci(:, ones(1, ProblemSize)) .* Ur*A';    %执行投影 
      yi = BoundConstraint_orgi(yi, Pop, LU);  %边界检查
      % yiFitness = feval(FHD, yi', varargin{1}); %评估适应度
      yiFitness = feval(fhd, yi', func);
      yiFitness = yiFitness'; %适应度
 %% To check stagnation
        % flag = false;
        bc_i_fit_var_old = bc_i_fit_var;
        for i = 1 : popsize
            NFES = NFES + 1;
            
            if yiFitness(i) < bc_i_fit_var
                bc_i_fit_var = yiFitness(i);
                bc_i_solution = yi(i, :);
                bc_i_index = i;
            end
            
            if NFES > Max_NFES
                break; 
            end
        end
        %需要把z清空，因为z是通过索引更新的
        % zi=[];

      % 
      % for i = 1 : popsize  %更换最好的适应度
      %     Nfes = Nfes + 1;
      %     if yiFitness(i) < BciFitVar
      %         BciFitVar = yiFitness(i);
      %         BciSolution = yi(i, :); 
      %     end
      % 
      %     if Nfes > MaxNfes; break; end
      % end 
      %%===================================================================
      dif = abs(fitness - yiFitness); %更新记忆
      I = (fitness > yiFitness); %存放优化成功的位置
      GdRD = rd(I == 1);   %好的rank
      GdCI = ci(I == 1);   %好的ci值
      DiffVal = dif(I == 1);   %成功个体的适应度差值

      Arch.pop=[Arch.pop;PopOld(I == 1, :)];
      Arch.OF=[Arch.OF;fitness(I==1)];
      nA=size(Arch.pop,1);
      % disp(nA==(Arch.NP+sum(I==1))); 


      % %% 自适应划分率
      % %比较两部分种群的成功率 这个有bug，会使得一部分种群直接为0
      % % disp(sum(I==1)==sum(I)); %测试
      % random_success_count = sum(I(ks));
      % best_success_count=sum(I(worst_index));
      % if n == 0
      %     random_successRate = 0;  % Or handle it differently as per the context
      % else
      %     random_successRate = random_success_count / n;
      % end
      % 
      % if (popsize - n) == 0
      %     best_successRate = 0;  % Or handle it differently as per the context
      % else
      %     best_successRate = best_success_count / (popsize - n);
      % end
      % 
      % % random_successRate=random_success_count/n;
      % 
      % % disp((popsize-n)==size(remaining_indices,1)); %测试
      % % best_successRate=best_success_count/(popsize-n);
      % if(random_successRate+best_successRate)==0
      %     M_p(M_iterator,1)=0.5;
      % else
      %     M_p(M_iterator,1)=random_successRate/(random_successRate+best_successRate);
      % end
      % 
      % M_iterator=M_iterator+1;
      % if(M_iterator>Msize)
      %     M_iterator=1;
      % end


      %%更新档案完毕
     

      %%
      %%===================================================================
      [fitness, I] = min([fitness, yiFitness], [], 2); %获得更好适应度的位置
      PopOld = Pop;
      PopOld(I == 2, :) = yi(I == 2, :);
      % PopOld2(I == 2, :) = Pop(I == 2, :); %保留被替换的父类个体

      






      %%===================================================================
      NumSucc = numel(GdRD); %统计成功的数量
      if NumSucc > 0 
	     SumDif = sum(DiffVal);
	     DiffVal = DiffVal / SumDif;
	     MemCI(MemI) =( (DiffVal' * (GdCI .^ 2)) / (DiffVal' * GdCI) +  MemCI(MemI))/2; 
         % MemCI(MemI) =( (DiffVal' * (GdCI .^ 2)) / (DiffVal' * GdCI) ); %使用这个
	     if max(GdRD) == 0 || MemRD(MemI)  == -1
	        MemRD(MemI)  = -1;
         else
	        MemRD(MemI) =( (DiffVal' * (GdRD .^ 2)) / (DiffVal' * GdRD) + MemRD(MemI) )/2;
            % MemRD(MemI) =( (DiffVal' * (GdRD .^ 2)) / (DiffVal' * GdRD)  ); %使用这个
         end
         % MemI = MemI + 1; %迭代器位置
	     % if MemI > Ms 
         %    MemI = 1; 
         % end
         MemI = mod(MemI, Ms-1) + 1;
      end

      %%===================================================================
      %%更新种群大小
      Plan_PopSize = round((((Min_Pop_Size - Max_Pop_Size) /Max_NFES) * NFES) + Max_Pop_Size);
      if popsize > Plan_PopSize
	     RedPop = popsize - Plan_PopSize;
	     if popsize - RedPop <  Min_Pop_Size
            RedPop = popsize - Min_Pop_Size;
         end
         popsize = popsize - RedPop;
	     for r = 1 : RedPop
	         [valBest indBest] = sort(fitness, 'ascend'); %%排序
	         worst_ind = indBest(end);
	         PopOld(worst_ind,:) = [];
             % PopOld2(worst_ind,:) = [];
	         Pop(worst_ind,:) = [];
	         fitness(worst_ind,:) = [];
             yiFitness(worst_ind,:)=[];
         end
	     Arch.NP = round(Ar * popsize); 
         % if size(Arch.pop, 1) > Arch.NP 
	     %    rndPos = randperm(size(Arch.pop, 1));
	     %    rndPos = rndPos(1 : Arch.NP);
	     %    Arch.pop = Arch.pop(rndPos, :);
	     % end
          % if size(Arch.pop, 1) > Arch.NP 
          
      end
      % 这个也应该每次都执行
      if nA > Arch.NP
          % Arch.pop=Arch.pop(:,(nA-Asize+1):nA);
          Arch.pop=Arch.pop((nA-Arch.NP+1):nA,:);
          %使用了ordering机制的话要更新适应度
          Arch.OF=Arch.OF((nA-Arch.NP+1):nA,:);
          % nA=Asize;
          % rndPos = randperm(size(Arch.pop, 1));
          % rndPos = rndPos(1 : Arch.NP);
          % Arch.pop = Arch.pop(rndPos, :);
      end


      % fitness_temp=fitness;
      

    
      % if rem(gg,10) == 1
      %     fprintf('best-so-far objective function at %d th iteration = %1.8e\n',gg,BciFitVar);
      % end
      % BciIndex(gg) = BciFitVar;



end  

      bc_i_error_val = bc_i_fit_var - Optm;
        if bc_i_error_val < val_2_reach
            bc_i_error_val = 0;
        end
        FF(run_id,func) = bc_i_error_val;
        fprintf('%d th run, best-so-far error value = %1.8e\n', run_id , bc_i_error_val)
        outcome = [outcome bc_i_error_val];  

       


       
        
       
        
        
        
       
    end %% end 1 run


    
end %% end 1 function run



end