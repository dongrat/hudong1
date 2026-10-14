function plot_sim_results(S, P)
%PLOT_SIM_RESULTS  Time / frequency plots for AK2401A two-point DMR sim
%  中文：绘制符号→RRC→两点路径→合成频偏与基带频谱。

t = S.t * 1e3;                    % ms
nShow = min(numel(t), round(0.04 * P.fs_Hz));  % ~40 ms window
idx = 1:nShow;

if P.saveFigs && ~exist(P.outDir, 'dir')
    mkdir(P.outDir);
end

%% Figure 1: symbols + ideal RRC deviation
figure('Name', 'DMR symbols and RRC', 'Color', 'w');
subplot(2,1,1);
stem((0:numel(S.symbols)-1)/P.Rs*1e3, S.symbols, 'filled', 'MarkerSize', 3);
grid on; ylabel('Symbol'); xlabel('t [ms]');
title('Random DMR 4FSK symbols \{-3,-1,+1,+3\}');
xlim([0, nShow/P.fs_Hz*1e3]);
subplot(2,1,2);
plot(t(idx), S.fIdeal(idx), 'k');
grid on; ylabel('Hz'); xlabel('t [ms]');
title(sprintf('RRC shaped deviation (alpha=%.2f, peak ~%.0f Hz)', ...
    P.alpha, max(abs(S.fIdeal))));
save_fig(P, '01_symbols_rrc');

%% Figure 2: Path A vs Path B commands and outputs
figure('Name', 'Two-point paths', 'Color', 'w');
subplot(3,1,1);
plot(t(idx), S.fA_cmd(idx), 'b', t(idx), S.fB_cmd(idx), 'r');
grid on; legend('Path A cmd (->OFST/SD)', 'Path B cmd (->DAC/VCO)', 'Location', 'northeast');
ylabel('Hz'); title('Two-point baseband / control commands');
subplot(3,1,2);
plot(t(idx), S.fPll(idx), 'b');
grid on; ylabel('Hz'); title(sprintf('Path A after PLL (f_{loop}=%.0f Hz)', P.fLoop_Hz));
subplot(3,1,3);
plot(t(idx), S.fVco(idx), 'r');
grid on; ylabel('Hz'); xlabel('t [ms]');
title(sprintf('Path B after DAC + SMF + Km (%.0f Hz/V)', P.Km_HzPerV));
save_fig(P, '02_twopoint_paths');

%% Figure 3: Combined vs ideal
figure('Name', 'Combined frequency trajectory', 'Color', 'w');
subplot(2,1,1);
plot(t(idx), S.fIdeal(idx), 'k', t(idx), S.fComb(idx), 'm', ...
     t(idx), S.fPll(idx), 'b--', t(idx), S.fVco(idx), 'r--');
grid on;
legend('Ideal RRC', 'Combined A+B', 'PLL only', 'VCO only', 'Location', 'northeast');
ylabel('Δf [Hz]'); title('Instantaneous frequency deviation');
subplot(2,1,2);
plot(t(idx), S.fComb(idx) - S.fIdeal(idx), 'Color', [0.6 0.2 0]);
grid on; ylabel('Hz'); xlabel('t [ms]');
title('Combined − Ideal (quantization / filter / ΣΔ residual)');
save_fig(P, '03_combined_vs_ideal');

%% Figure 4: Spectra
figure('Name', 'Spectra', 'Color', 'w');
[F, Pideal] = onesided_psd(S.fIdeal, P.fs_Hz);
[~, Pcomb]  = onesided_psd(S.fComb,  P.fs_Hz);
[~, Ppll]   = onesided_psd(S.fPll,   P.fs_Hz);
[~, Pvco]   = onesided_psd(S.fVco,   P.fs_Hz);
plot(F, 10*log10(Pideal+eps), 'k', ...
     F, 10*log10(Pcomb+eps),  'm', ...
     F, 10*log10(Ppll+eps),   'b', ...
     F, 10*log10(Pvco+eps),   'r');
grid on; xlim([0, min(8e3, P.fs_Hz/2)]);
xlabel('Frequency [Hz]'); ylabel('PSD [dB rel]');
legend('Ideal', 'Combined', 'Path A PLL', 'Path B VCO', 'Location', 'northeast');
title('Baseband frequency-deviation spectra');
save_fig(P, '04_spectra');

%% Figure 5: Baseband FM complex envelope spectrum (equivalent RF)
figure('Name', 'FM baseband spectrum', 'Color', 'w');
phi = 2*pi * cumsum(S.fComb) / P.fs_Hz;
s = exp(1j * phi);
[F2, Ps] = onesided_psd(s, P.fs_Hz);
plot(F2, 10*log10(Ps+eps), 'b');
grid on; xlim([0, min(12.5e3, P.fs_Hz/2)]);
xlabel('Frequency [Hz]'); ylabel('PSD [dB]');
title(sprintf('Equivalent FM complex envelope (|S(f)|), f_c = %.0f MHz implied', ...
    P.fRf_Hz/1e6));
save_fig(P, '05_fm_envelope_spectrum');

%% Figure 6: Architecture annotation (ASCII-only for Octave/gnuplot)
figure('Name', 'Architecture', 'Color', 'w');
axis off;
lines = {
    'AK2401A two-point TX (simulation)'
    sprintf('RRC alpha=%.2f | fLoop=%.0f Hz | Km=%.0f Hz/V | Kv=%.0f MHz/V', ...
        P.alpha, P.fLoop_Hz, P.Km_HzPerV, P.Kv_HzPerV/1e6)
    sprintf('Path A: OFST -> 18-bit SD Frac-N -> PLL LPF')
    sprintf('Path B: 12-bit DAC -> SMF(~20 kHz) -> VCO mod port')
    sprintf('Combined df = fPLL + fVCO | fPFD=%.1f MHz MOD=%d fs=%.0f kHz', ...
        P.fPfd_Hz/1e6, P.MOD, P.fs_Hz/1e3)
    sprintf('OFST LSB = fPFD/MOD = %.3f Hz', S.metaPll.ofstLsb_Hz)
    };
for k = 1:numel(lines)
    text(0.05, 0.9 - 0.12*(k-1), lines{k}, 'FontSize', 12, 'Interpreter', 'none');
end
save_fig(P, '06_architecture_notes');
end

function save_fig(P, name)
if ~P.saveFigs, return; end
fp = fullfile(P.outDir, [name '.png']);
try
    print(gcf, '-dpng', '-r120', fp);
catch
    % Octave/MATLAB print differences — ignore if headless fails
end
end

function [F, Pxx] = onesided_psd(x, fs)
x = x(:) - mean(x(:));
n = 2^nextpow2(min(numel(x), 8192));
X = fft(x(1:min(end,n)), n);
Pxx = abs(X(1:n/2+1)).^2 / (fs * n);
if n > 1
    Pxx(2:end-1) = 2*Pxx(2:end-1);
end
F = (0:n/2) * (fs/n);
end
