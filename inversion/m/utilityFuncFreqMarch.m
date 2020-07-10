function bispIndsFM = utilityFuncFreqMarch(L)
% Calculate the bispetral invariant indices that are non-zero
bispIndsFM = [];
for K=0:L
    c = 1;
    for l1=0:L
        for l2=0:l1
            for l=abs(l1-l2):min(l1+l2, L)
                if l1<=K && l2<=K && l<=K && (l1==K || l2==K || l==K)
                    bispIndsFM(K+1, c) = 1;
                    c = c + 1;
                    bispIndsFM(K+1, c) = 1;
                    c = c + 1;
                else
                    bispIndsFM(K+1, c) = 0;
                    c = c + 1;
                    bispIndsFM(K+1, c) = 0;
                    c = c + 1;
                end
            end
        end
    end
end
bispIndsFM = logical(bispIndsFM);
