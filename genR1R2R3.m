function [r1, r2, r3] = genR1R2R3(NP1, NP2, r0)
%要求r1 r2 r3 互不相等并且不能等于本身
NP0 = length(r0);
r3 = randperm(NP0);
for i = 1: 99999999
    Pos = (r3 == r0);
    if sum(Pos) == 0
        break;
    else
        r3 = randperm(NP0);
    end
     if i > 1000
        error('Can not genrate r3 in 1000 iterations');
     end
end
r1 = floor(rand(1, NP0) * NP1) + 1;

for i = 1 : 99999999
    Pos = (r1 == r0)|(r3 == r0);
    if sum(Pos) == 0
        break;
    else 
        r1(Pos) = floor(rand(1, sum(Pos)) * NP1) + 1;
    end
    if i > 1000
        error('Can not genrate r1 in 1000 iterations');
    end
end
r2 = floor(rand(1, NP0) * NP2) + 1;

for i = 1 : 99999999
    Pos = ((r2 == r1) | (r2 == r0)) | (r2 == r0);
    if sum(Pos)==0
        break;
    else 
        r2(Pos) = floor(rand(1, sum(Pos)) * NP2) + 1;
    end
    if i > 1000
        error('Can not genrate r2 in 1000 iterations');
    end
end
end