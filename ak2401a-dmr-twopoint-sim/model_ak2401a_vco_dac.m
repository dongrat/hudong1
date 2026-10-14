function [fVco_Hz, dacCodes, vSmf, meta] = model_ak2401a_vco_dac(fCmd_Hz, P)
%MODEL_AK2401A_VCO_DAC  Path B: freq cmd → 12-bit DAC → SMF → Kv → Hz
%
%  Datasheet: DAC+SMF (fc≈20 kHz) drives external VCO modulation pin (SMFOUT).
%  中文：路径 B — 12-bit DAC + 平滑滤波器驱动外部 VCO 调谐端。

fs = P.fs_Hz;
N = numel(fCmd_Hz);

% Map desired HF deviation to DAC voltage using modulation-port Km
% (not coarse Kv). v = f / Km ; full-scale uses typ Vpp
Km = P.Km_HzPerV;
vIdeal = (fCmd_Hz / Km) * P.pathBGain;
vPeak = P.dacVpp / 2;             % ±Vpp/2 around midscale
satMask = abs(vIdeal) > vPeak;
vClipped = max(min(vIdeal, vPeak), -vPeak);

% 12-bit mid-rise quantizer over ±vPeak
levels = 2^(P.dacBits - 1);
qStep = vPeak / levels;
dacCodes = round(vClipped / qStep);
dacCodes = max(min(dacCodes, levels-1), -levels);
vDac = dacCodes * qStep;

% SMF: approximate datasheet response (0 dB @1 kHz, −4 dB @20 kHz)
[b, a] = smf_lpf(P.smfFc_Hz, P.smfOrder, fs);
vSmf = filter(b, a, vDac);

% VCO modulation port: Δf = Km * v
fVco_Hz = vSmf * Km;

meta.qStep_V = qStep;
meta.b = b;
meta.a = a;
meta.vPeak = vPeak;
meta.satFraction = mean(satMask);
meta.Km_HzPerV = Km;
end

function [b, a] = smf_lpf(fc, order, fs)
% Simple Butterworth LPF; fc chosen so |H(j2π·20k)| ≈ −4 dB region
% Use 0.75*fc as -3dB to better match −4 dB @ 20 kHz for order-2
fc3 = 0.75 * fc;
Wn = min(fc3 / (fs/2), 0.99);
[b, a] = butter_local(order, Wn);
end

function [b, a] = butter_local(n, Wn)
% Minimal Butterworth lowpass (bilinear) for Octave/MATLAB without relying
% on signal-pkg edge cases. Supports n=1 or n=2.
if exist('butter', 'file') == 2
    [b, a] = butter(n, Wn);
    return;
end
% Fallback order-2 only
g = tan(pi*Wn/2);
if n == 1
    b = [g, g] / (1+g);
    a = [1, (g-1)/(g+1)];
else
    % Butterworth poles for n=2: angle 45°/135° → Q=1/sqrt(2)
    k = g*g;
    sqrt2 = sqrt(2);
    den = 1 + sqrt2*g + k;
    b = [k, 2*k, k] / den;
    a = [1, 2*(k-1)/den, (1 - sqrt2*g + k)/den];
end
end
