function shcOut = SHCProdFactor(shc, factor, L)
shcOut = zeros(size(shc));
for l=0:L
    for m=1:l
        shcOut(l*(l+1)+m+1) = shc(l*(l+1)+m+1)*factor;
        shcOut(l*(l+1)-m+1) = shc(l*(l+1)-m+1)*factor;
    end
end
end
