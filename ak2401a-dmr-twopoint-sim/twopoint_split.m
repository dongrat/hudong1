function [fA_cmd, fB_cmd, Hmeta] = twopoint_split(fIdeal_Hz, P)
%TWOPOINT_SPLIT  Complementary Path-A / Path-B baseband commands
%
%  Classic two-point:
%    Path A (PLL/ΣΔ): full deviation command → loop provides LPF
%    Path B (VCO DAC): complementary high-pass ≈ (1 − H_pll)
%  so H_A + H_B ≈ 1 across the RRC baseband.
%  中文：两点调制互补拆分 — A 路全带进 PLL，B 路高通进 VCO。

fs = P.fs_Hz;
wn = 2*pi*P.fLoop_Hz;
zeta = P.pllZeta;

% Same closed-loop LPF coefficients as PLL model (discrete)
[bL, aL] = pll_lpf_tf(wn, zeta, fs);

% DSP-side complementary split (common two-point firmware approach):
%   Path A gets LF content the PLL can track (into OFST / ΣΔ)
%   Path B gets HF residual (into DAC → VCO mod port)
fLp    = filter(bL, aL, fIdeal_Hz);
fA_cmd = fLp;
fB_cmd = fIdeal_Hz - fLp;

Hmeta.bL = bL;
Hmeta.aL = aL;
Hmeta.note = 'PathA_cmd=LPF(Ideal); PathB_cmd=Ideal-LPF(Ideal)';
end

function [b, a] = pll_lpf_tf(wn, zeta, fs)
T = 1/fs;
w = 2/T * tan(wn*T/2);
a2 = 1; a1 = 2*zeta*w; a0 = w^2;
A = a0*T^2 + 2*a1*T + 4*a2;
b = [a0*T^2, 2*a0*T^2, a0*T^2] / A;
a = [1, (2*a0*T^2 - 8*a2)/A, (a0*T^2 - 2*a1*T + 4*a2)/A];
end
