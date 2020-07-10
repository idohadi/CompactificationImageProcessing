function out = vec2structManopt(vecSHC, L)
out = struct();
for l=1:L
    out.(['M', num2str(l)]) = vecSHC(l*(l+1)+1 +(-l:l));
end
end
