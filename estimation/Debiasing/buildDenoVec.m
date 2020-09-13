function DenoVec = buildDenoVec(L)
% Builds the denoising vector
% 

DenoVec = []; % Denoising vector
c = 0;

for l1=0:L
    for l2=0:l1
        for l=abs(l1-l2):l1+l2
            c = c + 1;
            DenoVec(c) = 0;
            if l1==l && l2==0
                DenoVec(c) = DenoVec(c) + 2*l1 + 1;
            end
            if l2==l && l1==0
                DenoVec(c) = DenoVec(c) + 2*l2 + 1;
            end
        end
    end
end

DenoVec = DenoVec(:);
