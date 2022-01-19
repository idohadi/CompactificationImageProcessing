printBegEndMsg('Align inverted bispectrum and estimate SHC.', true);
[shcRelError, alignedSHCEst, optimalRotation, alignmentOutput] ...
    = alignSHC(shc, shcEstimator, bandlimit, tDesign);
save(fn, 'shcRelError', 'alignedSHCEst', 'optimalRotation', ...
    'alignmentOutput', '-append');
printBegEndMsg(num2str(shcRelError, ...
    'Align inverted bispectrum and estimate SHC.\n\tSHC rel err = %.3e.'), false);

printBegEndMsg('Inversion', false);


printBegEndMsg('Image estimation.', true);
imageEst ...
    = shc2image(alignedSHCEst, bandlimit, imageSize, ...
                interval, scalingParam);
imageEst = real(imageEst);
imageRelError = norm(im - imageEst, 'fro')/norm(im, 'fro');
imageBackProjection = shc2image(shc, bandlimit, imageSize, ...
                                interval, scalingParam);
imageRelErrorCleaned = imageRelError ...
    - norm(im - imageBackProjection)/norm(im, 'fro');
printBegEndMsg(num2str([imageRelError, imageRelErrorCleaned], ...
    'Image estimation.\n\tIm rel err = %.3e.\n\tIm rel err cln = %.3e.'), false);
save(fn, 'imageEst', 'imageRelError', 'imageRelErrorCleaned', '-append');
printBegEndMsg('Running test', true);

%% Produce figure
fig = figure;



savefig(fig, [fnNOEXT, '.fig']);

%% Shut down the diary
diary off;
