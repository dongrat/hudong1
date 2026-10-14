function h = rrc_filter(beta, span, sps)
%RRC_FILTER  Root-raised-cosine FIR impulse response (MATLAB/Octave)
%  h = rrc_filter(beta, span, sps)
%  beta: rolloff α, span: symbols, sps: samples/symbol
%  中文：根升余弦成形滤波器冲激响应（不依赖 rcosdesign）。

n = (span * sps) / 2;
t = (-n:n) / sps;                 % in symbol periods
h = zeros(size(t));

for k = 1:numel(t)
    ti = t(k);
    if abs(ti) < 1e-12
        h(k) = 1 - beta + 4*beta/pi;
    elseif beta > 0 && abs(abs(ti) - 1/(4*beta)) < 1e-12
        h(k) = (beta/sqrt(2)) * ( ...
            (1+2/pi)*sin(pi/(4*beta)) + (1-2/pi)*cos(pi/(4*beta)) );
    else
        num = sin(pi*ti*(1-beta)) + 4*beta*ti*cos(pi*ti*(1+beta));
        den = pi*ti*(1 - (4*beta*ti)^2);
        h(k) = num / den;
    end
end

h = h / sqrt(sum(h.^2));          % unit energy
end
