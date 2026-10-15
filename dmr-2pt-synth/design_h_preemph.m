function [b, a, info] = design_h_preemph(fs, smfB, smfA)
%DESIGN_H_PREEMPH  Band-limited inverse of the discrete SMF, for the H path.
%  The 3rd-order SMF is almost unity gain through 6 kHz and lags by about
%  15 deg at 3 kHz and 30 deg at 6 kHz. Its exact inverse has a triple pole
%  at z = -1 and runs away at Nyquist. This 2nd-order lead matches 1/Gsmf
%  from about 1 to 6 kHz (phase lead plus the small in-band boost) and
%  levels off near +14 dB, so the boost above the DMR band stays finite.
%
%  Analog prototype, DC gain 1, high-frequency gain (fp/fz)^2:
%    H(s) = (fp/fz)^2 * (s^2 + (wz/Qz) s + wz^2) / (s^2 + (wp/Qp) s + wp^2)
%  Corners were fitted at fs = 76800 Hz against the bilinear SMF. The loop
%  bandwidth and phase margin are not part of this fit.

if nargin < 1
    fs = 76800;
end

fz = 11952;
gainRatio = 2.220;     % wp/wz; high-frequency magnitude gain is gainRatio^2
Qz = 0.801;
Qp = 1.247;
wz = 2 * pi * fz;
wp = wz * gainRatio;
bb = (wp / wz)^2 * [1, wz / Qz, wz^2];
aa = [1, wp / Qp, wp^2];
[b, a] = tf_bilinear(bb, aa, fs);

info.fz_Hz = fz;
info.fp_Hz = fz * gainRatio;
info.Qz = Qz;
info.Qp = Qp;
info.hfGain_dB = 20 * log10(gainRatio^2);
info.b = b;
info.a = a;

probe = [3000, 6000];
for k = 1:numel(probe)
    [gdB, pd, pdSmf, gdProd, pdProd] = eval_at(b, a, smfB, smfA, fs, probe(k));
    info.(sprintf('f%d_Hz', probe(k))) = probe(k);
    info.(sprintf('boost%d_dB', probe(k) / 1000)) = gdB;
    info.(sprintf('lead%d_deg', probe(k) / 1000)) = pd;
    info.(sprintf('smf%d_deg', probe(k) / 1000)) = pdSmf;
    info.(sprintf('prod%d_dB', probe(k) / 1000)) = gdProd;
    info.(sprintf('prod%d_deg', probe(k) / 1000)) = pdProd;
end

fhi = (7000:200:fs / 2 - 50).';
zh = exp(1j * 2 * pi * fhi / fs);
Hh = polyval(b, zh) ./ polyval(a, zh);
[info.maxBoost_dB, ix] = max(20 * log10(abs(Hh)));
info.maxBoost_Hz = fhi(ix);
info.nyquist_dB = 20 * log10(abs(sum(b .* ((-1) .^ (0:numel(b) - 1))) ...
    / sum(a .* ((-1) .^ (0:numel(a) - 1)))));
end

function [pre_dB, pre_deg, smf_deg, prod_dB, prod_deg] = eval_at(b, a, smfB, smfA, fs, f0)
zz = exp(1j * 2 * pi * f0 / fs);
Hp = polyval(b, zz) / polyval(a, zz);
pre_dB = 20 * log10(abs(Hp));
pre_deg = angle(Hp) * 180 / pi;
if isempty(smfB)
    smf_deg = NaN;
    prod_dB = NaN;
    prod_deg = NaN;
    return;
end
Hs = polyval(smfB, zz) / polyval(smfA, zz);
Hc = Hp * Hs;
smf_deg = angle(Hs) * 180 / pi;
prod_dB = 20 * log10(abs(Hc));
prod_deg = angle(Hc) * 180 / pi;
end
