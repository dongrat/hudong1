function measure_fsk_acpr()
%MEASURE_FSK_ACPR  FSK error and 12.5 kHz ACPR of the existing two-point sim.
%  Does not change the modulator. Uses a longer capture than the demo.
%
%    cd /workspace/dmr-2pt-synth
%    measure_fsk_acpr

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
if exist('graphics_toolkit', 'file') || exist('graphics_toolkit', 'builtin')
    try
        graphics_toolkit('gnuplot');
    catch
    end
end

P = params_ak2404_dmr();
P.nSymbols = 400;   % a few hundred symbols; demo entry point stays at 96

fprintf('FSK / ACPR capture: %d symbols, seed %d\n', P.nSymbols, P.seed);
[~, symDev, mMod] = gen_dmr_4fsk(P);
nAcq = round(P.tAcquire_s * P.fsMod_Hz);
m = [zeros(1, nAcq), mMod];

tic;
fprintf('Simulating synthesizer (both paths, architecture unchanged)...\n');
out = sim_twopoint_pll(P, m, 'both');
fprintf('  %.1f s, DAC clips %d, 12-bit LSB %.4f Hz\n', toc, out.nDacClip, P.quantLsb_Hz);
mCmd = out.m_q;
dq = mCmd - m;
fprintf('  |m - m_q| max %.3f Hz (half of %.3f Hz LSB is %.3f Hz)\n', ...
    max(abs(dq)), P.quantLsb_Hz, P.quantLsb_Hz / 2);

sps = P.sps;
half = (P.rrcSpan / 2) * sps;
iValid0 = nAcq + half + 1;
iValid1 = nAcq + numel(mMod) - half;
% Pad so a ±1.5-symbol timing search stays inside real samples.
margin = 40;
i0 = iValid0 - margin;
i1 = iValid1 + margin;
inner = (iValid0:iValid1) - i0 + 1;

k0 = P.rrcSpan / 2 + 1;
k1 = P.nSymbols - P.rrcSpan / 2;
symIndex = k0:k1;
symSample = nAcq + (symIndex - 1) * sps + 1;
symLevel = symDev(symIndex);
symInSlice = symSample - i0 + 1;

fFloat = m(i0:i1);
fCmd = mCmd(i0:i1);
fSynth = out.fDevRf(i0:i1);

% Three fits. The reference is the second argument.
fprintf('Fitting synthesizer vs pre-quantization m...\n');
Rsm = fsk_fit(fSynth, fFloat, symInSlice, inner, P.fsMod_Hz);
fprintf('Fitting synthesizer vs m_q...\n');
Rsq = fsk_fit(fSynth, fCmd, symInSlice, inner, P.fsMod_Hz);
fprintf('Fitting m_q vs pre-quantization m...\n');
Rqm = fsk_fit(fCmd, fFloat, symInSlice, inner, P.fsMod_Hz);

fs = P.fsMod_Hz;
[acM, specM] = acpr_all(m(iValid0:iValid1), fs);
[acQ, specQ] = acpr_all(mCmd(iValid0:iValid1), fs);
[acS, specS] = acpr_all(out.fDevRf(iValid0:iValid1), fs);
acLeak = acpr_all(zeros(size(m(iValid0:iValid1))), fs);

print_fsk('synthesizer vs m(t)', Rsm);
print_fsk('synthesizer vs m_q', Rsq);
print_fsk('ideal FM of m_q vs m(t)', Rqm);
print_acpr('ideal FM of m(t)', acM);
print_acpr('ideal FM of m_q', acQ);
print_acpr('synthesizer', acS);
print_acpr('window leakage (unmodulated)', acLeak);

if ~exist(P.figDir, 'dir')
    mkdir(P.figDir);
end
if ~exist(P.mediaDir, 'dir')
    mkdir(P.mediaDir);
end
plot_fsk(P, Rsm, Rsq, Rqm, symLevel);
plot_acpr(P, specM, specQ, specS, acM, acQ, acS);
write_fsk_acpr(rootDir, P, Rsm, Rsq, Rqm, acM, acQ, acS, acLeak);
fprintf('Figures:\n  %s\n  %s\n', P.figDir, P.mediaDir);
end

function R = fsk_fit(fMeas, mRef, symInSlice, inner, fs)
% Minimize RMS of fMeas - f0 - k*mRef(t-tau) over tau, f0, k.
% tau > 0 delays the reference (mRef(t-tau)). Score only `inner` samples.
fMeas = fMeas(:);
mRef = mRef(:);
n = numel(fMeas);
span = 30;                 % samples; > 1 symbol at sps = 16, and < the 40-sample pad
coarse = -span:0.25:span;
[tauC, ~] = search_tau(fMeas, mRef, inner, coarse);
fine = (tauC - 0.5):0.02:(tauC + 0.5);
[tau, pack] = search_tau(fMeas, mRef, inner, fine);

mShift = delay_linear(mRef, tau);
keep = ~isnan(mShift);
f0 = pack.f0;
k = pack.k;
resid = nan(n, 1);
resid(keep) = fMeas(keep) - f0 - k * mShift(keep);
ref = k * mShift(keep);
err = resid(keep);
R.tau_s = tau / fs;
R.tau_samples = tau;
R.f0_Hz = f0;
R.k = k;
R.rmsErr_Hz = sqrt(mean(err.^2));
R.rmsRef_Hz = sqrt(mean(ref.^2));
R.fskWave_pct = 100 * R.rmsErr_Hz / R.rmsRef_Hz;
R.resid = resid;
R.t_s = (0:n - 1).' / fs;

% Symbol instants: measurement time t_sym+tau versus reference at t_sym.
tSym = symInSlice(:);
fAt = interp1((1:n).', fMeas, tSym + tau, 'linear', NaN);
mAt = mRef(tSym);
ok = ~isnan(fAt);
eSym = fAt(ok) - f0 - k * mAt(ok);
R.symErr_Hz = eSym;
R.fskSym_pct = 100 * sqrt(mean(eSym.^2)) / 1944;
R.symLevel = []; % filled by caller plot from the known levels, same ok mask
R.symOk = ok;
R.fAt = fAt;
R.mAt = mAt;
R.plot0 = inner(1);
end

function [tauBest, pack] = search_tau(fMeas, mRef, inner, taus)
best = inf;
tauBest = taus(1);
pack = struct('f0', 0, 'k', 1);
for i = 1:numel(taus)
    md = delay_linear(mRef, taus(i));
    keep = false(size(md));
    keep(inner) = ~isnan(md(inner));
    A = [ones(sum(keep), 1), md(keep)];
    coef = A \ fMeas(keep);
    err = fMeas(keep) - A * coef;
    e = sqrt(mean(err.^2));
    if e < best
        best = e;
        tauBest = taus(i);
        pack.f0 = coef(1);
        pack.k = coef(2);
    end
end
end

function y = delay_linear(x, tau)
n = numel(x);
y = interp1((1:n).', x(:), (1:n).' - tau, 'linear', NaN);
end

function [ac, spec] = acpr_all(fDev, fs)
fDev = fDev(:) - mean(fDev);
n = numel(fDev);
ph = 2 * pi * cumsum(fDev) / fs;
s = exp(1j * ph);
w = kaiser_local(n, 18);
X = fftshift(fft(s .* w));
pbin = abs(X).^2 / (n^2 * mean(w.^2));
f = ((-floor(n / 2):ceil(n / 2) - 1).') * (fs / n);
if numel(f) ~= numel(pbin)
    f = (-floor(n / 2):ceil(n / 2) - 1).' * (fs / n);
    f = f(1:numel(pbin));
end

Ptot = sum(pbin);
ac.Ptot = Ptot;
ac.window = window_ratios(f, pbin);
ac.mask = mask_ratios(f, pbin, Ptot);
ac.d2 = rect_ratios(f, pbin, Ptot, 4.25e3);

% Decimate the spectrum for plotting (keep peak in each bin group).
spec.f = f;
spec.p = pbin;
end

function r = window_ratios(f, pbin)
% Integrated power in ±6.25 kHz windows. Ratio to the on-carrier window.
half = 6.25e3;
centers = [-25e3, -12.5e3, 0, 12.5e3, 25e3];
names = {'m25', 'm12', 'c0', 'p12', 'p25'};
Pw = zeros(size(centers));
for i = 1:numel(centers)
    Pw(i) = sum(pbin(abs(f - centers(i)) <= half));
end
r.P = Pw;
r.dB = 10 * log10(Pw / Pw(3));
r.names = names;
end

function r = mask_ratios(f, pbin, Ptot)
% Nominal EN 300 113 annex C mask, centred on the adjacent/alternate channel.
% Power relative to total envelope power (constant-envelope carrier reference).
offsets = [-25e3, -12.5e3, 12.5e3, 25e3];
names = {'alt_m25', 'adj_m12', 'adj_p12', 'alt_p25'};
PdB = zeros(size(offsets));
for i = 1:numel(offsets)
    wgt = etsi_mask_power(f - offsets(i));
    PdB(i) = 10 * log10(sum(pbin .* wgt) / Ptot);
end
r.dB = PdB;
r.names = names;
r.adj_worst_dB = max(PdB(2:3));
r.alt_worst_dB = max(PdB([1, 4]));
end

function r = rect_ratios(f, pbin, Ptot, halfBw)
% Rectangular integration of width 2*halfBw, annex C.1.0 ±D2 alternative.
offsets = [-25e3, -12.5e3, 12.5e3, 25e3];
PdB = zeros(size(offsets));
for i = 1:numel(offsets)
    sel = abs(f - offsets(i)) <= halfBw;
    PdB(i) = 10 * log10(sum(pbin(sel)) / Ptot);
end
r.dB = PdB;
r.halfBw_Hz = halfBw;
end

function w = etsi_mask_power(df)
% Nominal 12.5 kHz mask. Attenuation breakpoints are distances from the
% measured-channel centre (EN 300 113 V3.1.1 annex C, tables C.2/C.3 and
% the 8.5 kHz necessary bandwidth in table 3).
x = abs(df) / 1e3;
bp = [0, 3.0, 4.25, 5.5, 9.5, 40];
at = [0, 2, 6, 26, 90, 90];
adb = interp1(bp, at, x, 'linear', 90);
w = 10 .^ (-adb / 10);
end

function w = kaiser_local(n, beta)
m = (0:n - 1).' - (n - 1) / 2;
t = 1 - (2 * m / (n - 1)).^2;
t(t < 0) = 0;
w = besseli(0, beta * sqrt(t)) / besseli(0, beta);
end

function print_fsk(label, R)
fprintf('%s\n', label);
fprintf('  waveform FSKerr %.3f %%   f0 %.3f Hz   k %.6f   tau %.3f us\n', ...
    R.fskWave_pct, R.f0_Hz, R.k, R.tau_s * 1e6);
fprintf('  symbol-instant FSKerr %.3f %%   (rms %.2f Hz / 1944 Hz)\n', ...
    R.fskSym_pct, sqrt(mean(R.symErr_Hz.^2)));
end

function print_acpr(label, ac)
fprintf('%s ACPR\n', label);
fprintf('  window ±6.25 kHz (dB vs on-channel window):');
fprintf(' %+.2f', ac.window.dB);
fprintf('\n');
fprintf('  EN mask vs total power (dBc) m25 m12 p12 p25:');
fprintf(' %.2f', ac.mask.dB);
fprintf('\n');
fprintf('    worst adjacent %.2f dBc, worst alternate %.2f dBc\n', ...
    ac.mask.adj_worst_dB, ac.mask.alt_worst_dB);
fprintf('  rect ±4.25 kHz vs total (dBc):');
fprintf(' %.2f', ac.d2.dB);
fprintf('\n');
end

function plot_fsk(P, Rsm, Rsq, Rqm, symLevel)
h = figure('visible', 'off');
set(h, 'position', [40 40 980 760]);
cM = [0.45 0.45 0.45];
cQ = [0.75 0.35 0.05];
cS = [0.12 0.25 0.65];

subplot(2, 1, 1);
i1 = min(numel(Rsm.resid), Rsm.plot0 + round(0.012 * P.fsMod_Hz) - 1);
idx = Rsm.plot0:i1;
tms = (idx - Rsm.plot0) / P.fsMod_Hz * 1e3;
plot(tms, Rqm.resid(idx), 'Color', cQ, 'LineWidth', 1.05);
hold on;
plot(tms, Rsq.resid(idx), 'Color', cS, 'LineWidth', 1.05);
plot(tms, Rsm.resid(idx), 'Color', cM, 'LineWidth', 1.0);
grid on;
xlabel('Time from start of fit window (ms)');
ylabel('Residual deviation (Hz)');
title(sprintf('Waveform fit   synth/m %.2f %%   synth/m_q %.2f %%   m_q/m %.2f %%', ...
    Rsm.fskWave_pct, Rsq.fskWave_pct, Rqm.fskWave_pct));
legend('m_q vs m(t)', 'synth vs m_q', 'synth vs m(t)', 'location', 'northeast');

subplot(2, 1, 2);
plot(symLevel(Rqm.symOk), Rqm.symErr_Hz, '.', 'Color', cQ, 'MarkerSize', 8);
hold on;
plot(symLevel(Rsq.symOk), Rsq.symErr_Hz, '.', 'Color', cS, 'MarkerSize', 8);
plot(symLevel(Rsm.symOk), Rsm.symErr_Hz, '.', 'Color', cM, 'MarkerSize', 7);
plot([-2200 2200], [0 0], 'k-', 'LineWidth', 0.7);
grid on;
xlabel('Known symbol deviation (Hz)');
ylabel('Symbol-centre residual (Hz)');
title(sprintf('Symbol centres / 1944 Hz   synth/m %.2f %%   synth/m_q %.2f %%   m_q/m %.2f %%', ...
    Rsm.fskSym_pct, Rsq.fskSym_pct, Rqm.fskSym_pct));
legend('m_q vs m(t)', 'synth vs m_q', 'synth vs m(t)', 'location', 'northeast');
save_png(h, P, '06_fsk_error.png');
end

function plot_acpr(P, specM, specQ, specS, acM, acQ, acS)
h = figure('visible', 'off');
set(h, 'position', [40 40 980 640]);
ref = max(specM.p);
dbM = 10 * log10(specM.p / ref + 1e-30);
dbQ = 10 * log10(specQ.p / ref + 1e-30);
dbS = 10 * log10(specS.p / ref + 1e-30);
plot(specM.f / 1e3, dbM, 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0);
hold on;
plot(specQ.f / 1e3, dbQ, 'Color', [0.75 0.35 0.05], 'LineWidth', 1.0);
plot(specS.f / 1e3, dbS, 'Color', [0.1 0.25 0.65], 'LineWidth', 1.0);
for fk = [-25, -12.5, 12.5, 25]
    plot([fk fk], [-120 5], '--', 'Color', [0.75 0.2 0.15], 'LineWidth', 0.8);
end
for fk = [-31.25, -18.75, -6.25, 6.25, 18.75, 31.25]
    plot([fk fk], [-120 5], ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.6);
end
grid on;
xlim([-40 40]);
ylim([-110 5]);
xlabel('Offset from 450 MHz carrier (kHz)');
ylabel('Normalized envelope power (dB)');
title(sprintf('Mask adjacent   m %.1f   m_q %.1f   synth %.1f dBc', ...
    acM.mask.adj_worst_dB, acQ.mask.adj_worst_dB, acS.mask.adj_worst_dB));
legend('ideal FM of m(t)', 'ideal FM of m_q', 'synthesizer', 'location', 'northeast');
save_png(h, P, '07_acpr_spectrum.png');
end

function save_png(h, P, name)
out = fullfile(P.figDir, name);
print(h, out, '-dpng', '-r140');
close(h);
copyfile(out, fullfile(P.mediaDir, name));
end

function write_fsk_acpr(rootDir, P, Rsm, Rsq, Rqm, acM, acQ, acS, acLeak)
fid = fopen(fullfile(rootDir, 'fsk_acpr_results.txt'), 'w');
fprintf(fid, 'nSymbols %d\n', P.nSymbols);
fprintf(fid, 'seed %d\n', P.seed);
fprintf(fid, 'sps %d\n', P.sps);
fprintf(fid, 'fsMod_Hz %.0f\n', P.fsMod_Hz);
fprintf(fid, 'fPfd_Hz %.0f\n', P.fPfd_Hz);
fprintf(fid, 'INT %d\n', P.INT);
fprintf(fid, 'FRAC %d\n', P.FRAC);
fprintf(fid, 'MOD %d\n', P.MOD);
fprintf(fid, 'quant_bits %d\n', P.quantBits);
fprintf(fid, 'quant_lsb_Hz %.8f\n', P.quantLsb_Hz);
write_fit(fid, 'synth_vs_m', Rsm);
write_fit(fid, 'synth_vs_mq', Rsq);
write_fit(fid, 'mq_vs_m', Rqm);
write_ac(fid, 'm', acM);
write_ac(fid, 'mq', acQ);
write_ac(fid, 'synth', acS);
fprintf(fid, 'leak_window_dB');
fprintf(fid, ' %.4f', acLeak.window.dB);
fprintf(fid, '\n');
fprintf(fid, 'leak_mask_dB');
fprintf(fid, ' %.4f', acLeak.mask.dB);
fprintf(fid, '\n');
fclose(fid);
end

function write_fit(fid, tag, R)
fprintf(fid, '%s_fsk_wave_pct %.6f\n', tag, R.fskWave_pct);
fprintf(fid, '%s_f0_Hz %.6f\n', tag, R.f0_Hz);
fprintf(fid, '%s_k %.8f\n', tag, R.k);
fprintf(fid, '%s_tau_us %.6f\n', tag, R.tau_s * 1e6);
fprintf(fid, '%s_fsk_sym_pct %.6f\n', tag, R.fskSym_pct);
fprintf(fid, '%s_sym_rms_Hz %.6f\n', tag, sqrt(mean(R.symErr_Hz.^2)));
end

function write_ac(fid, tag, ac)
fprintf(fid, '%s_window_dB', tag);
fprintf(fid, ' %.4f', ac.window.dB);
fprintf(fid, '\n');
fprintf(fid, '%s_mask_dB', tag);
fprintf(fid, ' %.4f', ac.mask.dB);
fprintf(fid, '\n');
fprintf(fid, '%s_rect_dB', tag);
fprintf(fid, ' %.4f', ac.d2.dB);
fprintf(fid, '\n');
end
