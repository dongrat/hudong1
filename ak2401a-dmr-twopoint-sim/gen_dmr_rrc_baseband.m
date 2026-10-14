function [t, fIdeal_Hz, symbols, info] = gen_dmr_rrc_baseband(P)
%GEN_DMR_RRC_BASEBAND  Random DMR 4FSK symbols → RRC → freq deviation [Hz]
%  中文：随机 DMR 4FSK 符号经 RRC 成形，得到理想瞬时频偏轨迹。

if isfield(P, 'seed') && ~isempty(P.seed)
    rand('seed', P.seed); %#ok<RAND>
end

symbols = P.symLevels(randi(numel(P.symLevels), 1, P.numSymbols));

% Upsample to sps pulses then RRC filter
ups = zeros(1, P.numSymbols * P.sps);
ups(1:P.sps:end) = symbols;

h = rrc_filter(P.alpha, P.rrcSpan, P.sps);
% Normalize so constant symbol=+1 → mean amplitude 1 after shaping
% (then ×devPerLevel maps steady ±3 → ±1944 Hz). RRC alone has ISI, so
% instantaneous peaks may slightly exceed ±1944 Hz — that is expected.
probe = zeros(1, 64*P.sps); probe(1:P.sps:end) = 1;
yprobe = conv(probe, h, 'same');
ssGain = mean(yprobe(16*P.sps : 48*P.sps));
h = h / ssGain;

shaped = conv(ups, h, 'same');

% Map |±3| → ±1944 Hz, |±1| → ±648 Hz (linear in symbol amplitude)
fIdeal_Hz = shaped * P.devPerLevel;
t = (0:numel(fIdeal_Hz)-1) / P.fs_Hz;

info.hRrc = h;
info.peakDev_Hz = max(abs(fIdeal_Hz));
info.Rs = P.Rs;
info.fs = P.fs_Hz;
info.sps = P.sps;
end
