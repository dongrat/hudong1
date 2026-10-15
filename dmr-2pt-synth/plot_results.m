function plot_results(P, R)
%PLOT_RESULTS  Write the five required PNG figures.

figDir = P.figDir;
mediaDir = P.mediaDir;
if ~exist(figDir, 'dir')
    mkdir(figDir);
end
if ~exist(mediaDir, 'dir')
    mkdir(mediaDir);
end

try
    graphics_toolkit('gnuplot');
catch
end

tms = R.t * 1e3;
i0 = R.iEx0;
i1 = R.iEx1;
tex = tms(i0:i1);

% ---- 1. baseband m(t) and symbol levels ----
h = figure('visible', 'off');
set(h, 'position', [40 40 980 520]);
plot(tex, R.m(i0:i1), 'Color', [0.1 0.1 0.1], 'LineWidth', 1.3);
hold on;
symT = R.tSym_ms;
symY = R.symDev;
% Keep symbols that fall inside the excerpt (time base is absolute).
sel = symT >= tex(1) - 1e-6 & symT <= tex(end) + 1e-6;
stem(symT(sel), symY(sel), 'filled', 'Color', [0.15 0.4 0.75], 'MarkerSize', 5);
lev = [1944, 648, -648, -1944];
for k = 1:4
    plot([tex(1) tex(end)], [lev(k) lev(k)], '--', 'Color', [0.75 0.75 0.75]);
end
grid on;
xlabel('Time (ms)');
ylabel('Frequency deviation (Hz)');
title('DMR 4FSK symbols and SRRC deviation m(t)');
legend('m(t) SRRC', 'dibit deviation', 'location', 'northeast');
ylim([-3200 3200]);
save_png(h, figDir, mediaDir, '01_baseband_mt_symbols.png');

% ---- 2. each path's contribution ----
h = figure('visible', 'off');
set(h, 'position', [40 40 980 560]);
plot(tex, R.m(i0:i1), 'Color', [0.45 0.45 0.45], 'LineWidth', 1.4);
hold on;
plot(tex, R.pll.fDevRf(i0:i1), 'Color', [0.8 0.2 0.15], 'LineWidth', 1.2);
plot(tex, R.vco.fDevRf(i0:i1), 'Color', [0.1 0.45 0.75], 'LineWidth', 1.2);
plot(tex, R.pll.fDevRf(i0:i1) + R.vco.fDevRf(i0:i1), 'Color', [0.15 0.55 0.25], 'LineWidth', 1.05);
grid on;
xlabel('Time (ms)');
ylabel('RF frequency deviation (Hz)');
title('Path contributions across the ~1 kHz PLL bandwidth');
legend('ideal m(t)', 'PLL path only (lowpass)', 'VCO path only (highpass)', ...
    'sum of the two paths', 'location', 'northeast');
save_png(h, figDir, mediaDir, '02_path_contributions.png');

% ---- 3. combined vs ideal, with lock ----
h = figure('visible', 'off');
set(h, 'position', [40 40 980 720]);
subplot(2, 1, 1);
acq = 1:R.nAcq;
plot(tms(acq), R.both.fDevRf(acq) / 1e3, 'Color', [0.1 0.2 0.55], 'LineWidth', 1.2);
grid on;
xlabel('Time (ms)');
ylabel('RF deviation (kHz)');
title('Acquisition: initial VCO error, then lock (modulation off)');
subplot(2, 1, 2);
plot(tex, R.m(i0:i1), 'Color', [0.55 0.55 0.55], 'LineWidth', 1.6);
hold on;
plot(tex, R.both.fDevRf(i0:i1), 'Color', [0.05 0.05 0.05], 'LineWidth', 1.15);
grid on;
xlabel('Time (ms)');
ylabel('RF frequency deviation (Hz)');
title('Combined instantaneous deviation vs ideal m(t)');
legend('ideal m(t)', 'synthesizer output', 'location', 'northeast');
save_png(h, figDir, mediaDir, '03_combined_vs_ideal.png');

% ---- 4. deviation error ----
h = figure('visible', 'off');
set(h, 'position', [40 40 980 520]);
errEx = R.both.fDevRf(i0:i1) - R.m(i0:i1);
plot(tex, errEx, 'Color', [0.55 0.1 0.45], 'LineWidth', 1.1);
grid on;
xlabel('Time (ms)');
ylabel('Deviation error (Hz)');
title(sprintf('Deviation error    RMS %.1f Hz over the settled symbols (peak %.1f Hz)', ...
    R.rmsErr_Hz, R.peakErr_Hz));
save_png(h, figDir, mediaDir, '04_deviation_error.png');

% ---- 5. output spectrum ----
h = figure('visible', 'off');
set(h, 'position', [40 40 980 720]);
iv0 = R.iValid0;
iv1 = R.iValid1;
fs = P.fsMod_Hz;
fDev = R.both.fDevRf(iv0:iv1);
mDev = R.m(iv0:iv1);
[fAxis, pSim] = envelope_spectrum(fDev, fs);
[~, pId] = envelope_spectrum(mDev, fs);
subplot(2, 1, 1);
plot(fAxis / 1e3, 10 * log10(pId / max(pId) + 1e-18), 'Color', [0.55 0.55 0.55], 'LineWidth', 1.3);
hold on;
plot(fAxis / 1e3, 10 * log10(pSim / max(pId) + 1e-18), 'Color', [0.1 0.2 0.55], 'LineWidth', 1.15);
grid on;
xlim([-20 20]);
ylim([-90 5]);
xlabel('Offset from 450 MHz carrier (kHz)');
ylabel('Normalized power (dB)');
title('TX complex-envelope spectrum');
legend('ideal m(t)', 'synthesizer', 'location', 'northeast');

[fD, pMd] = dev_spectrum(mDev, fs);
[~, pSd] = dev_spectrum(fDev, fs);
subplot(2, 1, 2);
plot(fD / 1e3, 10 * log10(pMd / max(pMd) + 1e-18), 'Color', [0.55 0.55 0.55], 'LineWidth', 1.3);
hold on;
plot(fD / 1e3, 10 * log10(pSd / max(pMd) + 1e-18), 'Color', [0.8 0.25 0.15], 'LineWidth', 1.1);
grid on;
xlim([-20 20]);
ylim([-70 5]);
xlabel('Frequency (kHz)');
ylabel('Normalized PSD (dB)');
title('Spectrum of the frequency-deviation waveform');
legend('ideal m(t)', 'synthesizer deviation', 'location', 'northeast');
save_png(h, figDir, mediaDir, '05_output_spectrum.png');
end

function save_png(h, figDir, mediaDir, name)
out = fullfile(figDir, name);
set(h, 'paperpositionmode', 'auto');
print(h, out, '-dpng', '-r140');
close(h);
copyfile(out, fullfile(mediaDir, name));
end

function [f, p] = envelope_spectrum(fDev, fs)
fDev = fDev(:);
n = numel(fDev);
ph = 2 * pi * cumsum(fDev) / fs;
s = exp(1j * ph);
w = 0.5 - 0.5 * cos(2 * pi * (0:n - 1).' / (n - 1));
S = fftshift(fft(s .* w));
p = abs(S).^2;
p = p / sum(w.^2);
f = ((-n / 2):(n / 2 - 1)).' * (fs / n);
% n even: fftshift bin 0 is at index n/2+1, and the negative end is one longer
% if n is even, the vector from -n/2 to n/2-1 has length n. Good.
if numel(f) ~= numel(p)
    f = (-floor(n / 2):ceil(n / 2) - 1).' * (fs / n);
end
end

function [f, p] = dev_spectrum(x, fs)
x = x(:) - mean(x);
n = numel(x);
w = 0.5 - 0.5 * cos(2 * pi * (0:n - 1).' / (n - 1));
S = fftshift(fft(x .* w));
p = abs(S).^2;
p = p / sum(w.^2);
f = (-floor(n / 2):ceil(n / 2) - 1).' * (fs / n);
end
