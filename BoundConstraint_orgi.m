function zi = BoundConstraint_orgi (zi, Pop, LU)
[NP, ~] = size(Pop); 
xl = repmat(LU(1, :), NP, 1);
Pos = zi < xl;
zi(Pos) = (Pop(Pos) + xl(Pos)) / 2;
xu = repmat(LU(2, :), NP, 1);
Pos = zi > xu;
zi(Pos) = (Pop(Pos) + xu(Pos)) / 2;
end