function out = struct2vecManopt(structSHC, L)
out = zeros((L+1)^2, 1);
out(1) = 1;
for l=1:L
    out(l*(l+1)+1 +(-l:l)) = structSHC.(['M', num2str(l)]);
end
end
