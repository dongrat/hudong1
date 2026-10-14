%RUN_AK2401A_DMR_TWOPOINT  Main entry — AK2401A DMR two-point TX simulation
%
%  Architecture (datasheet-aligned):
%    Random DMR 4FSK → RRC (α=0.2) → two-point split
%      Path A → Frequency Offset (OFST) → 18-bit ΣΔ Frac-N PLL
%      Path B → 12-bit DAC → SMF (~20 kHz) → external VCO (Kv)
%    Combined instantaneous deviation + equivalent FM baseband spectrum
%
%  Usage (MATLAB or Octave):
%    cd <this folder>
%    run_ak2401a_dmr_twopoint
%
%  中文：AK2401A 两点调制 + DMR 4FSK 频综仿真入口脚本。

clear; close all; clc;

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);

% Octave: load signal package if present (butter)
if exist('OCTAVE_VERSION', 'builtin')
    try
        pkg load signal;
    catch
        warning('octave-signal not loaded; using local LPF fallbacks where needed');
    end
    % Prefer visible figures when display exists; otherwise save PNG
    try, graphics_toolkit('gnuplot'); catch, end
end

%% 1) Parameters
P = params_ak2401a_dmr();
if P.verbose
    fprintf('=== AK2401A DMR two-point TX simulation ===\n');
    fprintf('Rs=%.0f baud, α=%.2f, sps=%d, fs=%.1f kHz\n', ...
        P.Rs, P.alpha, P.sps, P.fs_Hz/1e3);
    fprintf('fRF=%.1f MHz, fPFD=%.2f MHz, INT=%d, FRAC=%d, MOD=%d\n', ...
        P.fRf_Hz/1e6, P.fPfd_Hz/1e6, P.INT, P.FRAC, P.MOD);
    fprintf('ICP≈%.1f uA, fLoop=%.0f Hz, Kv=%.1f MHz/V, Km=%.0f Hz/V, DAC %d-bit @ %.0f kHz\n', ...
        P.Icp_uA, P.fLoop_Hz, P.Kv_HzPerV/1e6, P.Km_HzPerV, P.dacBits, P.dacFs_Hz/1e3);
end

%% 2) Random DMR symbols → RRC → ideal Δf [Hz]
[t, fIdeal, symbols, infoSrc] = gen_dmr_rrc_baseband(P);

%% 3) Two-point complementary split
[fA_cmd, fB_cmd, Hmeta] = twopoint_split(fIdeal, P);

%% 4) Path A: OFST / ΣΔ / PLL
[fPll, ofstCodes, metaPll] = model_ak2401a_sd_pll(fA_cmd, P);

%% 5) Path B: DAC / SMF / VCO
[fVco, dacCodes, vSmf, metaVco] = model_ak2401a_vco_dac(fB_cmd, P);

%% 6) Combined RF frequency trajectory (baseband FM equivalent)
fComb = fPll + fVco;
% Ignore filter startup (~2 symbol periods) for error stats
i0 = max(1, 2*P.sps);
errRms = sqrt(mean((fComb(i0:end) - fIdeal(i0:end)).^2));
errPeak = max(abs(fComb(i0:end) - fIdeal(i0:end)));

if P.verbose
    fprintf('OFST LSB = %.4f Hz\n', metaPll.ofstLsb_Hz);
    fprintf('Peak |ideal|=%.1f Hz, combined RMS error=%.2f Hz, peak err=%.1f Hz\n', ...
        infoSrc.peakDev_Hz, errRms, errPeak);
    fprintf('OFST code range: [%d, %d]\n', min(ofstCodes), max(ofstCodes));
    fprintf('DAC code range:  [%d, %d]  (sat fraction=%.2f%%)\n', ...
        min(dacCodes), max(dacCodes), 100*metaVco.satFraction);
end

%% 7) Pack + plot
S = struct();
S.t = t;
S.symbols = symbols;
S.fIdeal = fIdeal;
S.fA_cmd = fA_cmd;
S.fB_cmd = fB_cmd;
S.fPll = fPll;
S.fVco = fVco;
S.fComb = fComb;
S.ofstCodes = ofstCodes;
S.dacCodes = dacCodes;
S.vSmf = vSmf;
S.metaPll = metaPll;
S.metaVco = metaVco;
S.Hmeta = Hmeta;
S.errRms = errRms;

plot_sim_results(S, P);

%% 8) Save workspace snapshot for post-analysis
outMat = fullfile(thisDir, 'last_run_results.mat');
try
    save(outMat, 'P', 'S', 'symbols', 'fIdeal', 'fPll', 'fVco', 'fComb', 't');
    if P.verbose
        fprintf('Saved %s\n', outMat);
        if P.saveFigs
            fprintf('Figures → %s\n', P.outDir);
        end
    end
catch ME
    warning('Could not save MAT file: %s', ME.message);
end

fprintf('Done.\n');
