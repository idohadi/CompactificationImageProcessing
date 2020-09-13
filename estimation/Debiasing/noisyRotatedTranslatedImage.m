function im_noisy = noisyRotatedTranslatedImage(im, N, sigma)

bin_bnd = 30;

[h, w] = size(im);

im_noisy = cell(N, 1);
for j=1:N
    sample = imrotate(im, randi(360), 'bicubic', 'crop');
    
    trans_directions = binornd(2*bin_bnd, 0.5, [1, 2]) - bin_bnd;
    sample = imtranslate(sample, trans_directions, 'cubic', 'OutputView', 'same');
    
    im_noisy{j} = sample + sqrt(sigma)*randn(h, w);
end
