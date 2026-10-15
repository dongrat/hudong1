function P = params_ak2404_dmr()
%PARAMS_AK2404_DMR  AK2404 two-point TX synthesizer + DMR 4FSK parameters.
%  Numbers taken from AK2404 datasheet 220900015-E-04 are marked "DS".
%  Anything the datasheet does not specify is marked "ASSUMED".

%% DMR 4FSK (ETSI TS 102 361-1) — specified by the task, not the AK2404 DS
P.Rs = 4800;                  % symbols/s
P.alpha = 0.2;                % SRRC roll-off
P.rrcSpan = 12;               % symbols (truncation; 12 keeps constant-symbol ripple small)
P.nSymbols = 96;
P.seed = 2404;
% Dibit map: 01 +1944, 00 +648, 10 -648, 11 -1944 Hz
P.dev_Hz = struct('d01', 1944, 'd00', 648, 'd10', -648, 'd11', -1944);

%% Reference, PFD, fractional-N
% 19.6608 MHz is not one of the three recommended TCXO frequencies
% (18.432 / 19.2 / 24.576). Section 13.6: TCXO_SET = 0 or 3 turns the
% CLOCK RATE CONVERTER off and uses the TCXO frequency as-is, so
% fMCLK = fTCXO. 19.6608 is nearer 18.432, so the bypass setting is
% TCXO_SET = 3. fMCLK is then 19.6608 MHz, not the table's nominal 18.432.
% Datasheet note 47: a master clock other than 18.432 or 24.576 changes
% digital timing. fPFD stays inside the 1–25 MHz range. PLL_R = 1.
P.fTcxo_Hz = 19.6608e6;
P.fMclk_Hz = P.fTcxo_Hz;      % converter bypassed
P.tcxoSet = 3;                % converter off (closer of the two bypass codes)
P.pllR = 1;
P.fPfd_Hz = P.fTcxo_Hz / P.pllR;
P.fRf_Hz = 450e6;             % DS BAND1 test frequency; DMR UHF 12.5 kHz grid
P.divSel = 2;                 % DS BAND1 local divider (DIVSEL = 2)
P.denRatio = 2;               % DS: offset step ratio tracks the LO divider
P.fVco_Hz = P.fRf_Hz * P.divSel;
% Exact carrier: 19.6608e6 * (45 + 795/1024) = 900e6.
% The 125 Hz example raster (MOD = fPFD/125) does not hit 900 MHz exactly.
% The DSM modulus is the largest multiple of 1024 that still fits in 23 bits
% (MOD <= 8388607). FRAC scales by the same factor, so the ratio is unchanged.
% MOD = 1024 itself would quantize OFS_MDLT onto fPFD/1024 = 19.2 kHz at the
% VCO (9.6 kHz at RF) and would put a fractional spur on the alternate channel.
P.INT = 45;
P.modReduced = 1024;
P.fracReduced = 795;
P.modScale = floor(8388607 / P.modReduced);   % 8191
P.MOD = P.modReduced * P.modScale;            % 8387584
P.FRAC = P.fracReduced * P.modScale;          % 6511845
P.ofsDen = 8388607;           % DS max unsigned 23-bit OFS_DEN
P.ofsBits = 16;               % DS OFS_MDLT signed width
P.quantBits = 12;             % only the upper 12 bits are used
P.quantShift = 16;            % one 12-bit count = 16 OFS LSBs
P.ofsLsb_Hz = P.fPfd_Hz / P.ofsDen;
P.quantLsb_Hz = P.quantShift * P.ofsLsb_Hz;   % ~37.5 Hz at RF
P.dsmOrder = 3;               % ASSUMED MASH 1-1-1 (DS says 23-bit ΔΣ, not the order)
P.dsmAtInt = 1;               % DS: normally 1 when offset modulation is used

%% Charge pump (DS formula; CPFINE register is ASSUMED)
P.Rbias2_kohm = 27;           % DS: connect 27 kΩ to BIAS2
P.CPFINE = 7;                 % ASSUMED (0..31). Eval-board value is not given.
P.Icp_uA = (2160 / P.Rbias2_kohm) * (P.CPFINE + 1);
P.Icp_A = P.Icp_uA * 1e-6;

%% VCO and loop (BW target is DS; Kvco and the R/C realization are ASSUMED)
P.Kvco_HzPerV = 10e6;         % ASSUMED external-VCO tune gain
P.fLoop_Hz = 1000;            % DS eval board: loop BW approximately 1 kHz
P.phaseMargin_deg = 55;       % ASSUMED
P.fPole3_Hz = 200e3;          % ASSUMED 3rd pole (DS: 3rd-order LF, frac noise > 1 MHz)
P.fDevInit_Hz = 400e3;        % ASSUMED initial VCO frequency error (shows lock)

%% MODDAC / SMF / VCO modulation port (DS electricals; Km ASSUMED)
P.sps = 16;                   % SRRC samples per symbol
P.fsMod_Hz = P.Rs * P.sps;    % 16 * 4800 = 76800 Hz
P.modFs_Hz = P.fsMod_Hz;      % MODFS; DS maximum is 96 kHz
P.dacBits = 12;               % DS
P.dacVpp = 1.35;              % DS typ full-scale at MODDAC_AG = 00
P.dacVref = 1.45;             % DS typ reference level (bias, not deviation)
P.modDacAg_dB = 0;            % DS setting 00
P.modDacDg = 511;             % DS digital gain 0 dB
P.Km_HzPerV = 20e3;           % ASSUMED VCO modulation-port sensitivity (Hz at VCO)
P.smfFcAnalog_Hz = [];        % filled below
P.smfOrder = 3;               % ASSUMED; fits DS typ -4 dB @ 20 kHz and -44 dB @ 100 kHz

%% Simulation schedule
P.tAcquire_s = 0.006;
P.nPfdPerMod = P.fPfd_Hz / P.fsMod_Hz;   % 19.6608e6 / 76800 = 256 exactly

if abs(P.nPfdPerMod - 256) > 1e-9
    error('expected 256 PFD cycles per modulation sample, got %.6f', P.nPfdPerMod);
end
P.nPfdPerMod = 256;
fCar = P.fPfd_Hz * (P.INT + P.FRAC / P.MOD);
if abs(fCar - P.fVco_Hz) > 1e-3
    error('carrier is not exact: VCO %.6f Hz', fCar);
end

%% Derived loop filter and SMF
L = design_loop_filter(P);
P.loop = L;
[P.smfB, P.smfA, P.smfFcAnalog_Hz] = smf_discrete(P.fsMod_Hz);

%% Paths
P.figDir = fullfile(fileparts(mfilename('fullpath')), 'figures');
P.mediaDir = '/cursor/stores/self/media';
end

function [b, a, fc] = smf_discrete(fs)
% 3rd-order Butterworth whose analog gain is -4 dB at 20 kHz (DS typ).
ratio = (10^(4 / 10) - 1)^(1 / 6);
fc = 20e3 / ratio;
wp = 2 * fs * tan(pi * fc / fs);   % prewarp the cutoff
% (s^2 + s + 1)*(s + 1) = s^3 + 2 s^2 + 2 s + 1, scaled by wp
b = wp^3;
a = [1, 2 * wp, 2 * wp^2, wp^3];
[b, a] = tf_bilinear(b, a, fs);
end
