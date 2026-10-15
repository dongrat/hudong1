function [bits, symDev, m, dibit] = gen_dmr_4fsk(P)
%GEN_DMR_4FSK  Fixed-seed random dibits -> SRRC frequency deviation m(t).
%  m is sampled at P.fsMod_Hz. DC gain is unity: a constant dibit settles
%  at its symbol deviation.

rng(P.seed, 'twister');
bits = randi([0, 1], 1, P.nSymbols * 2);

symDev = zeros(1, P.nSymbols);
dibit = zeros(1, P.nSymbols);
for k = 1:P.nSymbols
    b1 = bits(2 * k - 1);
    b0 = bits(2 * k);
    dibit(k) = b1 * 2 + b0;
    if b1 == 0 && b0 == 1
        symDev(k) = P.dev_Hz.d01;
    elseif b1 == 0 && b0 == 0
        symDev(k) = P.dev_Hz.d00;
    elseif b1 == 1 && b0 == 0
        symDev(k) = P.dev_Hz.d10;
    else
        symDev(k) = P.dev_Hz.d11;
    end
end

h = srrc_filter(P.alpha, P.rrcSpan, P.sps);
ups = zeros(1, P.nSymbols * P.sps);
ups(1:P.sps:end) = symDev;
m = conv(ups, h, 'same');
end
