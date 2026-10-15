function run_dmr_2pt_synth()
%RUN_DMR_2PT_SYNTH  AK2404 two-point DMR synthesizer simulation.
%  Octave or MATLAB, from this directory:
%    run_dmr_2pt_synth
%  Headless Octave:
%    octave --no-gui --eval "cd /workspace/dmr-2pt-synth; run_dmr_2pt_synth"

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
if exist('graphics_toolkit', 'builtin') || exist('graphics_toolkit', 'file')
    try
        graphics_toolkit('gnuplot');
    catch
    end
end

P = params_ak2404_dmr();
check_bilinear();
check_srrc_dc(P);

fprintf('AK2404 two-point DMR synthesizer\n');
fprintf('  RF %.3f MHz   VCO %.3f MHz   DIV %d   N = %.6f (INT %d FRAC %d / MOD %d)\n', ...
    P.fRf_Hz / 1e6, P.fVco_Hz / 1e6, P.divSel, P.loop.N, P.INT, P.FRAC, P.MOD);
fprintf('  fPFD %.3f MHz   ICP %.1f uA (CPFINE %d, BIAS2 %g kOhm)\n', ...
    P.fPfd_Hz / 1e6, P.Icp_uA, P.CPFINE, P.Rbias2_kohm);
fprintf('  Loop %.0f Hz, PM %.0f deg, C1 %.3f nF, C2 %.3f nF, R2 %.1f ohm, pole3 %.0f kHz\n', ...
    P.fLoop_Hz, P.phaseMargin_deg, P.loop.C1_F * 1e9, P.loop.C2_F * 1e9, ...
    P.loop.R2_ohm, P.fPole3_Hz / 1e3);
fprintf('  Kvco %.2f MHz/V   Km %.1f kHz/V   SMF fc %.2f kHz (analog fit)\n', ...
    P.Kvco_HzPerV / 1e6, P.Km_HzPerV / 1e3, P.smfFcAnalog_Hz / 1e3);
print_smf_gain(P);

[bits, symDev, mMod, dibit] = gen_dmr_4fsk(P);
nAcq = round(P.tAcquire_s * P.fsMod_Hz);
m = [zeros(1, nAcq), mMod];
t = (0:numel(m) - 1) / P.fsMod_Hz;

fprintf('  Symbols %d, seed %d, acquire %.1f ms, MOD ticks %d x %d PFD cycles\n', ...
    P.nSymbols, P.seed, P.tAcquire_s * 1e3, numel(m), P.nPfdPerMod);
fprintf('  First dibits:');
for k = 1:8
    fprintf(' %d%d', bits(2 * k - 1), bits(2 * k));
end
fprintf('\n');

tic;
fprintf('Simulating both paths...\n');
R.both = sim_twopoint_pll(P, m, 'both');
fprintf('  %.1f s   divider word %d..%d   DAC clips %d\n', toc, ...
    R.both.nInstMin, R.both.nInstMax, R.both.nDacClip);
tic;
fprintf('Simulating PLL path only...\n');
R.pll = sim_twopoint_pll(P, m, 'pll');
fprintf('  %.1f s\n', toc);
tic;
fprintf('Simulating VCO path only...\n');
R.vco = sim_twopoint_pll(P, m, 'vco');
fprintf('  %.1f s\n', toc);

% Plots and the error numbers use m_q, the shared 16-bit command.
R.m = R.both.m_q;
R.mFloat = m;
R.t = t;
R.symDev = symDev;
R.bits = bits;
R.dibit = dibit;
R.nAcq = nAcq;
R.mMod = mMod;

sps = P.sps;
half = (P.rrcSpan / 2) * sps;
R.iValid0 = nAcq + half + 1;
R.iValid1 = nAcq + numel(mMod) - half;
nExSym = 12;
R.iEx0 = R.iValid0;
R.iEx1 = min(R.iValid1, R.iValid0 + nExSym * sps - 1);

% Symbol impulse times on the absolute simulation clock.
R.tSym_ms = (P.tAcquire_s + (0:P.nSymbols - 1) / P.Rs) * 1e3;

iv = R.iValid0:R.iValid1;
err = R.both.fDevRf(iv) - R.m(iv);
R.rmsErr_Hz = sqrt(mean(err.^2));
R.peakErr_Hz = max(abs(err));
sumErr = (R.pll.fDevRf(iv) + R.vco.fDevRf(iv)) - R.both.fDevRf(iv);
R.rmsLin_Hz = sqrt(mean(sumErr.^2));
R.lock_s = find_lock(R.both.fDevRf(1:nAcq), P.fsMod_Hz);

% Settled acquisition residual (last 1 ms of the lock interval).
nTail = round(0.001 * P.fsMod_Hz);
R.acqResidual_Hz = sqrt(mean(R.both.fDevRf(nAcq - nTail + 1:nAcq).^2));

fprintf('Lock time: ');
if isnan(R.lock_s)
    fprintf('not detected inside the acquisition window\n');
else
    fprintf('%.2f ms\n', R.lock_s * 1e3);
end
fprintf('Acquisition residual RMS (last 1 ms): %.2f Hz\n', R.acqResidual_Hz);
fprintf('Deviation error vs m_q: RMS %.2f Hz, peak %.2f Hz\n', R.rmsErr_Hz, R.peakErr_Hz);
fprintf('Quantizer |m - m_q| max %.3f Hz (LSB %.4f Hz)\n', ...
    max(abs(R.mFloat - R.m)), P.quantLsb_Hz);
fprintf('Linearity (pll-only + vco-only - both) RMS: %.2f Hz\n', R.rmsLin_Hz);

if ~exist(P.figDir, 'dir')
    mkdir(P.figDir);
end
plot_results(P, R);
write_metrics(rootDir, P, R);

fprintf('Figures written to\n  %s\n  %s\n', P.figDir, P.mediaDir);
end

function check_bilinear()
fs = 96e3;
fc = 1e3;
wc = 2 * pi * fc;
[b, a] = tf_bilinear(wc, [1, wc], fs);
z = exp(1j * 2 * pi * fc / fs);
H = polyval(fliplr(b), z) / polyval(fliplr(a), z);
gdb = 20 * log10(abs(H));
if abs(gdb + 3) > 0.15
    error('bilinear self-check failed: gain at fc is %.2f dB', gdb);
end
end

function check_srrc_dc(P)
h = srrc_filter(P.alpha, P.rrcSpan, P.sps);
nSym = 40;
ups = zeros(1, nSym * P.sps);
ups(1:P.sps:end) = 1944;
m = conv(ups, h, 'same');
mid = m(round(end * 0.4):round(end * 0.6));
fprintf('SRRC DC check: constant +1944 mean %.2f Hz, ripple %.1f Hz pp\n', ...
    mean(mid), max(mid) - min(mid));
if abs(mean(mid) - 1944) > 1
    error('SRRC normalization error');
end
end

function print_smf_gain(P)
for f = [1e3, 20e3]
    z = exp(1j * 2 * pi * f / P.fsMod_Hz);
    H = polyval(fliplr(P.smfB), z) / polyval(fliplr(P.smfA), z);
    fprintf('  SMF discrete gain at %g kHz: %.2f dB\n', f / 1e3, 20 * log10(abs(H)));
end
end

function tLock = find_lock(dev, fs)
holdN = round(0.0005 * fs);
thresh = 2000;
tLock = NaN;
for i = 1:(numel(dev) - holdN)
    if all(abs(dev(i:i + holdN - 1)) < thresh)
        tLock = (i - 1) / fs;
        return;
    end
end
end

function write_metrics(rootDir, P, R)
fid = fopen(fullfile(rootDir, 'results_summary.txt'), 'w');
fprintf(fid, 'carrier_Hz %.0f\n', P.fRf_Hz);
fprintf(fid, 'vco_Hz %.0f\n', P.fVco_Hz);
fprintf(fid, 'fPfd_Hz %.0f\n', P.fPfd_Hz);
fprintf(fid, 'INT %d\n', P.INT);
fprintf(fid, 'FRAC %d\n', P.FRAC);
fprintf(fid, 'MOD %d\n', P.MOD);
fprintf(fid, 'Icp_uA %.3f\n', P.Icp_uA);
fprintf(fid, 'C1_nF %.6f\n', P.loop.C1_F * 1e9);
fprintf(fid, 'C2_nF %.6f\n', P.loop.C2_F * 1e9);
fprintf(fid, 'R2_ohm %.6f\n', P.loop.R2_ohm);
fprintf(fid, 'R3_ohm %.6f\n', P.loop.R3_ohm);
fprintf(fid, 'C3_nF %.6f\n', P.loop.C3_F * 1e9);
fprintf(fid, 'Kvco %.6e\n', P.Kvco_HzPerV);
fprintf(fid, 'Km %.6e\n', P.Km_HzPerV);
fprintf(fid, 'smf_fc %.6f\n', P.smfFcAnalog_Hz);
fprintf(fid, 'lock_s %.8f\n', R.lock_s);
fprintf(fid, 'acq_rms_Hz %.6f\n', R.acqResidual_Hz);
fprintf(fid, 'rms_err_Hz %.6f\n', R.rmsErr_Hz);
fprintf(fid, 'peak_err_Hz %.6f\n', R.peakErr_Hz);
fprintf(fid, 'lin_rms_Hz %.6f\n', R.rmsLin_Hz);
fprintf(fid, 'nInst %d %d\n', R.both.nInstMin, R.both.nInstMax);
fprintf(fid, 'dac_clip %d\n', R.both.nDacClip);
fprintf(fid, 'seed %d\n', P.seed);
fclose(fid);
end
