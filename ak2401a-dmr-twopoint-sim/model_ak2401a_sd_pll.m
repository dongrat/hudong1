function [fPll_Hz, ofstCodes, meta] = model_ak2401a_sd_pll(fCmd_Hz, P)
%MODEL_AK2401A_SD_PLL  Path A: OFST → 18-bit ΣΔ fractional-N → PLL LPF
%
%  Models AK2401A digital FM via Frequency Offset (§13.7.3):
%    Δf = fPFD * OFST / MOD
%  then closed-loop PLL tracking (type-II approx as 2nd-order LPF at fLoop).
%  Optional coarse ΣΔ quantization noise on the fractional command.
%  中文：路径 A — 频偏寄存器进 ΣΔ/分频器，再经 PLL 闭环低通。

fs = P.fs_Hz;
N = numel(fCmd_Hz);

% --- Quantize command to OFST codes (18-bit 2's complement) ---
ofstReal = fCmd_Hz * (P.MOD / P.fPfd_Hz);
ofstCodes = round(ofstReal);
ofstMax = 2^(18-1) - 1;           % signed 18-bit
ofstCodes = max(min(ofstCodes, ofstMax), -ofstMax-1);
fQuant_Hz = ofstCodes * (P.fPfd_Hz / P.MOD);

% --- Lightweight ΣΔ: 2nd-order error feedback on fractional residue ---
% Emulates that only integer N-divider pulses leave the DSM each PFD cycle.
% At baseband fs << fPFD we inject shaped quantization in the baseband model.
fracErr = ofstReal - ofstCodes;
e1 = 0; e2 = 0;
sdNoise = zeros(1, N);
for n = 1:N
    u = fracErr(n) + 2*e1 - e2;
    q = round(u);
    e2 = e1;
    e1 = u - q;
    sdNoise(n) = q - fracErr(n);
end
% Scale SD noise to Hz and attenuate (fs model, not full-rate PFD)
fSd_Hz = fQuant_Hz + sdNoise * (P.fPfd_Hz / P.MOD) * 0.05;

% --- Closed-loop PLL dynamics ---
wn = 2*pi*P.fLoop_Hz;
zeta = P.pllZeta;
if isfield(P, 'pllPreSplit') && P.pllPreSplit
    % Command is already LPF'd to ~fLoop: model residual tracking lag with a
    % faster 1st-order (≈5×fLoop) so we do not double-apply the loop LPF.
    tau = 1 / (2*pi*5*P.fLoop_Hz);
    a1 = exp(-1/(fs*tau));
    b = [1-a1];
    a = [1, -a1];
    fPll_Hz = filter(b, a, fSd_Hz);
    meta.mode = 'presplit-light-track';
else
    % Full command into ΣΔ: closed-loop provides the two-point LPF
    [b, a] = butter_from_2nd_order(wn, zeta, fs);
    fPll_Hz = filter(b, a, fSd_Hz);
    meta.mode = 'full-command-lpf';
end

meta.b = b;
meta.a = a;
meta.fQuant_Hz = fQuant_Hz;
meta.ofstLsb_Hz = P.fPfd_Hz / P.MOD;
meta.formula = 'f = fPFD*(INT + (FRAC+OFST)/MOD)';
end

function [b, a] = butter_from_2nd_order(wn, zeta, fs)
% Bilinear transform of continuous  ωn^2 / (s^2 + 2ζωn s + ωn^2)
T = 1/fs;
% Prewarp
w = 2/T * tan(wn*T/2);
a2 = 1;
a1 = 2*zeta*w;
a0 = w^2;
% bilinear s = 2/T *(1-z^-1)/(1+z^-1)
A = a0*T^2 + 2*a1*T + 4*a2;
b0 = a0*T^2 / A;
b1 = 2*a0*T^2 / A;
b2 = a0*T^2 / A;
aa0 = 1;
aa1 = (2*a0*T^2 - 8*a2) / A;
aa2 = (a0*T^2 - 2*a1*T + 4*a2) / A;
b = [b0 b1 b2];
a = [aa0 aa1 aa2];
end
