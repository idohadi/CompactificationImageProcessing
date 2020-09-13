function im_noisy = noisyImages(im, N, sigma)


[im_h, im_w] = size(im);
im_noisy = cell(N, 1);
for n=1:N
    im_noisy{n} = im + sqrt(sigma)*randn(im_h, im_w);
end
