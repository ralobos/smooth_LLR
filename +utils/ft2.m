function out = ft2(x)

out = fftshift(fft2(ifftshift(x)));