function [R] = localrelief(bed,nc)

[ny,nx] = size(bed);
%nc = 4;
%nc = 2;
% resRelief ~ 2nc*res
R = zeros(ny,nx);

for i = (1+nc):(ny-nc),
    
    for j = (1+nc):(nx-nc),
        
        M = bed((i-nc):(i+nc),(j-nc):(j+nc));
        R(i,j) = max(M(:)) - min(M(:));
        
    end;
    
end;

