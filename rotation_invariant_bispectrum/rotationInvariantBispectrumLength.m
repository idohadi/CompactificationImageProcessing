function bLen = rotationInvariantBispectrumLength(truncation, angularLimits)
%%
% TODO docs

%% Compute the length of the bispectrum vector
angularLimits = angularLimits(2:end);

bLen = 0;
for k1=1:truncation
    for k2=1:min(k1, truncation-k1)
        bLen = bLen ...
            + angularLimits(k1) ...
                * angularLimits(k2) ...
                * angularLimits(k1+k2);
    end
end
