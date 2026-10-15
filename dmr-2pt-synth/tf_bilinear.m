function [bz, az] = tf_bilinear(b, a, fs)
%TF_BILINEAR  Bilinear transform (s = 2*fs*(1-z^-1)/(1+z^-1)).
%  b, a are continuous-time coefficients, highest power of s first.
%  bz, az are discrete coefficients for filter(), b0 + b1*z^-1 + ...

b = b(:).';
a = a(:).';
M = max(numel(a), numel(b)) - 1;
if numel(a) < M + 1
    a = [zeros(1, M + 1 - numel(a)), a];
end
if numel(b) < M + 1
    b = [zeros(1, M + 1 - numel(b)), b];
end

bz = zeros(1, M + 1);
az = zeros(1, M + 1);
for k = 0:M
    sk = M - k;
    pk = 1;
    for i = 1:sk
        pk = conv(pk, [1, -1]);
    end
    for i = 1:(M - sk)
        pk = conv(pk, [1, 1]);
    end
    scale = (2 * fs)^sk;
    bz = bz + b(k + 1) * scale * pk;
    az = az + a(k + 1) * scale * pk;
end
bz = bz / az(1);
az = az / az(1);
end
