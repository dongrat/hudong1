function h = srrc_filter(beta, span, sps)
%SRRC_FILTER  Square-root raised-cosine impulse response.
%  span is in symbols. The taps are normalized so a constant symbol
%  stream of amplitude A (impulses every sps samples) settles to A.

n = -span * sps / 2 : span * sps / 2;
t = n / sps;
h = zeros(size(t));
for i = 1:numel(t)
    ti = t(i);
    fourBt = 4 * beta * ti;
    if abs(ti) < 1e-12
        h(i) = 1 - beta + 4 * beta / pi;
    elseif abs(abs(fourBt) - 1) < 1e-8
        h(i) = (beta / sqrt(2)) * ((1 + 2 / pi) * sin(pi / (4 * beta)) ...
            + (1 - 2 / pi) * cos(pi / (4 * beta)));
    else
        h(i) = (sin(pi * ti * (1 - beta)) + 4 * beta * ti * cos(pi * ti * (1 + beta))) ...
            / (pi * ti * (1 - fourBt^2));
    end
end
h = h / sum(h) * sps;
end
