function [r1, r2,r3]=gnR1R2_prob2(r0,prob,Indexes)
%% 请注意，用此函数，种群是没有进行排序的
% gnA1A2 generate two column vectors r1 and r2 of size NP1 & NP2, respectively
%    r1's elements are choosen from {1, 2, ..., NP1} & r1(i) ~= r0(i)
%    r2's elements are choosen from {1, 2, ..., NP2} & r2(i) ~= r1(i) & r2(i) ~= r0(i)
%
% Call:
%    [r1 r2 ...] = gnA1A2(NP1)   % r0 is set to be (1:NP1)'
%    [r1 r2 ...] = gnA1A2(NP1, r0) % r0 should be of length NP1
%
% Version: 2.1  Date: 2008/07/01
% Written by Jingqiao Zhang (jingqiao@gmail.com)

NP0 = length(r0); %种群大小
%% %% 生成r1，r1是压力选择的随机数
%r1 = floor(rand(1, NP0) * NP1) + 1;  
temp=randsample(r0,NP0,true,prob);
r1=Indexes(temp); %%实际所得到的索引
%for i = 1 : inf
for i = 1 : 99999999
    pos = (r1 == r0);  %检查有没有选到自己的情况
    if sum(pos) == 0  %如果没有就退出
        break;
    else % regenerate r1 if it is equal to r0 %有就重新生成
        %r1(pos) = floor(rand(1, sum(pos)) * NP1) + 1;
        temp=randsample(r0,sum(pos),true,prob);
        r1(pos)=Indexes(temp);
        %r1(pos)=randsample(r0,sum(pos),true,prob);
    end
    if i > 1000 % this has never happened so far
        error('Can not genrate r1 in 1000 iterations');
    end
end
% if(archiveSize==0)
%     r2=[];
% else
%     r2 = floor(rand(1, NP0) * archiveSize) + 1;
%     for i = 1 : inf
%     for i = 1 : 99999999
%         pos = ((r2 == r1) | (r2 == r0));
%         if sum(pos)==0
%             break;
%         else % regenerate r2 if it is equal to r0 or r1 同样要保证不能选到自身且不能等于r1
%             r2(pos) = floor(rand(1, sum(pos)) * archiveSize) + 1;
%         end
%         if i > 1000 % this has never happened so far
%             error('Can not genrate r2 in 1000 iterations');
%         end
%     end
% end



temp=randsample(r0,NP0,true,prob);
r2=Indexes(temp);
%for i = 1 : inf
for i = 1 : 99999999
    pos = ((r2 == r1) | (r2 == r0));
    if sum(pos)==0
        break;
    else % regenerate r2 if it is equal to r0 or r1 同样要保证不能选到自身且不能等于r1
        temp=randsample(r0,sum(pos),true,prob);
        r2(pos) =Indexes(temp);
    end
    if i > 1000 % this has never happened so far
        error('Can not genrate r2 in 1000 iterations');
    end
end

r3 = floor(rand(1, NP0) * NP0) + 1;

for i = 1 : 99999999
    Pos = (r3 == r0)|(r3 == r1)|(r3==r2);
    if sum(Pos) == 0
        break;
    else 
        r3(Pos) = floor(rand(1, sum(Pos)) * NP0) + 1;
    end
    if i > 1000
        error('Can not genrate r3 in 1000 iterations');
    end
end

end




