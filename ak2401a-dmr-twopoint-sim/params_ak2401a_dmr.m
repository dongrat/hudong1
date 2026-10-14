function P = params_ak2401a_dmr()
%PARAMS_AK2401A_DMR  AK2401A + DMR two-point TX simulation parameters
% 中文：集中存放数据手册摘录参数与 DMR / 仿真假设。

%% ---- DMR physical layer (ETSI TS 102 361) ----
P.Rs          = 4800;           % symbol rate [baud]
P.bitsPerSym  = 2;              % 4FSK
P.alpha       = 0.2;            % RRC rolloff
P.rrcSpan     = 8;              % filter span [symbols]
P.symLevels   = [-3 -1 1 3];    % 4FSK symbol alphabet
P.devPerLevel = 648;            % Hz per |symbol step| → ±648 / ±1944 Hz
P.numSymbols  = 256;            % random source length
P.seed        = 20260314;

%% ---- AK2401A reference / fractional-N (datasheet) ----
P.fTcxo_Hz    = 19.2e6;         % recommended TCXO
P.fPfd_Hz     = 2.4e6;          % PFD (datasheet OFST examples use 2.4 MHz)
P.MOD         = 262143;         % prefer max 18-bit MOD
P.fRf_Hz      = 450e6;          % UHF DMR carrier (external VCO / LO)
% Integer/frac for carrier (before modulation offset)
Ntot          = P.fRf_Hz / P.fPfd_Hz;
P.INT         = floor(Ntot);
P.FRAC        = round((Ntot - P.INT) * P.MOD);
if P.FRAC >= P.MOD
    P.FRAC = P.MOD - 1;
end

%% ---- Charge pump (datasheet formulas) ----
P.Rbias2_kOhm = 27;             % typ BIAS2
P.CPFINE      = 7;              % register 0..31 → ICP = ICP_MIN*(n+1)
P.Icp_uA      = (2160 / P.Rbias2_kOhm) * (P.CPFINE + 1);

%% ---- External VCO + loop (NOT in datasheet — assumptions) ----
P.Kv_HzPerV   = 20e6;           % coarse VCO tune gain [Hz/V] (PLL loop design)
% Modulation-port sensitivity for SMFOUT → VCO mod pin (much smaller than Kv):
% sized so ~±2 kHz uses a healthy fraction of 1.35 Vpp DAC swing.
P.Km_HzPerV   = 8e3;            % [Hz/V]  headroom for HF residual vs 1.35 Vpp
P.fLoop_Hz    = 800;            % closed-loop BW for Path-A LPF model
P.pllZeta     = 0.707;          % damping of 2nd-order closed-loop approx

%% ---- Path B: 12-bit DAC + SMF (datasheet) ----
P.dacBits     = 12;
P.dacFs_Hz    = 96e3;           % <= 200 kHz max
P.dacVpp      = 1.35;           % typ SMFOUT swing [Vpp]
P.smfFc_Hz    = 20e3;           % SMF −4 dB @ 20 kHz
P.smfOrder    = 2;              % approx LPF order for SMF

%% ---- Simulation rates ----
P.sps         = round(P.dacFs_Hz / P.Rs);  % samples per symbol @ DAC rate
P.fs_Hz       = P.Rs * P.sps;             % baseband / control sample rate
assert(P.fs_Hz <= 200e3, 'fs exceeds AK2401A DAC max 200 kHz');

%% ---- Two-point split ----
% Path A: PLL tracks low-freq OFST; Path B: complementary high-pass via DAC→VCO
P.pathBGain   = 1.0;            % VCO-path gain trim (1.0 = ideal match)
% DSP already applies LPF for Path A command; PLL model uses light tracking
% (near-unity inside band) instead of applying fLoop LPF a second time.
P.pllPreSplit = true;

%% ---- Plots / misc ----
P.saveFigs    = true;
P.outDir      = fullfile(fileparts(mfilename('fullpath')), 'figures');
P.verbose     = true;

end
