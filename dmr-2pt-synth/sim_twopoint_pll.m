function out = sim_twopoint_pll(P, mRf, mode)
%SIM_TWOPOINT_PLL  Behavioral AK2404 two-point fractional-N synthesizer.
%  One PFD period per iteration: MASH divider, tri-state PFD, charge pump,
%  3rd-order loop filter, VCO, plus the MODDAC/SMF direct-modulation path.
%
%  mode 'both' : Lch OFS_MDLT and Hch MODOUT both carry m(t)
%       'pll'  : low-frequency path only (closed-loop lowpass)
%       'vco'  : direct VCO path only (closed-loop highpass)
%
%  mRf is the desired RF frequency deviation in Hz at fsMod (76.8 kHz).
%  It is quantized to the upper 12 bits of the OFS word before either path:
%    q16 = round(m * OFS_DEN / fPFD)
%    q12 = clamp(round(q16 / 16), -2048, 2047)
%    m_q = q12 * 16 * fPFD / OFS_DEN
%  Lch OFS_MDLT = q12 * 16. Hch is rounded from m_q only:
%    Hch = clamp(round(m_q * DIVSEL / Km * 4096 / 1.35), -2048, 2047)
%  out.fDevRf is the TX-divider output deviation, (fVCO - fVCO0)/DIV.

if nargin < 3
    mode = 'both';
end

usePll = strcmp(mode, 'both') || strcmp(mode, 'pll');
useVco = strcmp(mode, 'both') || strcmp(mode, 'vco');

nTick = numel(mRf);
nPer = P.nPfdPerMod;
fPfd = P.fPfd_Hz;
T = 1 / fPfd;
Kv = P.Kvco_HzPerV;
Icp = P.Icp_A;
Ctot = P.loop.Ctot_F;
T2 = P.loop.T2;
e1 = exp(-T / P.loop.T1);
e3 = exp(-T / P.loop.T3);
INT = P.INT;
MOD = P.MOD;
fFree = P.fVco_Hz;

% Both paths see the same quantized deviation. Do this before either channel.
[m_q, qWord, q12] = quantize_ofs(P, mRf);

% --- Hch: 12-bit MODDAC from m_q only, then SMF and Km (Hz at the VCO) ---
fDirTick = zeros(1, nTick);
nClip = 0;
if useVco
    % code = m_q * DIVSEL / Km * 4096 / Vpp, signed 12-bit.
    codeUnc = round(m_q * P.divSel / P.Km_HzPerV * 4096 / P.dacVpp);
    code = min(2047, max(-2048, codeUnc));
    nClip = sum(code ~= codeUnc);
    vAc = code * (P.dacVpp / 4096);
    vSmf = filter(P.smfB, P.smfA, vAc(:));
    fDirTick = (P.Km_HzPerV * vSmf.') * 10^(P.modDacAg_dB / 20);
end

% --- Lch: OFS_MDLT = q12*16, then the fraction added to FRAC/MOD ---
accWord = zeros(1, nTick);
nCarry = zeros(1, nTick);
if usePll
    ofs = qWord;
    % delta = DEN_RATIO * OFS / OFS_DEN ; VCO frequency step = fPFD * delta
    add = round(P.denRatio * ofs * MOD / P.ofsDen);
    acc = P.FRAC + add;
    carry = zeros(1, nTick);
    hi = acc >= MOD;
    while any(hi)
        acc(hi) = acc(hi) - MOD;
        carry(hi) = carry(hi) + 1;
        hi = acc >= MOD;
    end
    lo = acc < 0;
    while any(lo)
        acc(lo) = acc(lo) + MOD;
        carry(lo) = carry(lo) - 1;
        lo = acc < 0;
    end
    accWord = acc;
    nCarry = carry;
else
    accWord(:) = P.FRAC;
end

% Initial tune voltage = initial frequency error / Kvco (filter at rest).
v = P.fDevInit_Hz / Kv;
q = v * Ctot;
p = 0;
w = 0;
pe = 0;

acc1 = 0;
acc2 = 0;
acc3 = 0;
c2z = 0;
c3z = 0;
c3zz = 0;

fMean = zeros(1, nTick);
nMin = INT;
nMax = INT;

for it = 1:nTick
    ain = accWord(it);
    nBase = INT + nCarry(it);
    fd = fDirTick(it);
    accF = 0;
    for j = 1:nPer
        acc1 = acc1 + ain;
        c1 = 0;
        if acc1 >= MOD
            acc1 = acc1 - MOD;
            c1 = 1;
        end
        acc2 = acc2 + acc1;
        c2 = 0;
        if acc2 >= MOD
            acc2 = acc2 - MOD;
            c2 = 1;
        end
        acc3 = acc3 + acc2;
        c3 = 0;
        if acc3 >= MOD
            acc3 = acc3 - MOD;
            c3 = 1;
        end
        mash = c1 + (c2 - c2z) + (c3 - 2 * c3z + c3zz);
        c3zz = c3z;
        c3z = c3;
        c2z = c2;

        ninst = nBase + mash;
        if ninst < 4
            ninst = 4;
        end
        if ninst < nMin
            nMin = ninst;
        elseif ninst > nMax
            nMax = ninst;
        end

        fv = fFree + Kv * v + fd;
        % Time error of one divider edge: Tdiv - Tref, in reference cycles.
        % Locks to fVCO = fPFD * mean(n_inst) (arithmetic mean of the MASH word).
        pe = pe + (ninst * fPfd / fv - 1);
        if pe > 1
            pe = 1;
        elseif pe < -1
            pe = -1;
        end

        icp = Icp * pe;
        w = e3 * w + (1 - e3) * icp;
        p = e1 * p + (1 - e1) * w;
        q = q + T * p;
        v = (q + T2 * p) / Ctot;
        fv = fFree + Kv * v + fd;
        accF = accF + fv;
    end
    fMean(it) = accF / nPer;
end

out.fVco = fMean;
out.fDevRf = (fMean - P.fVco_Hz) / P.divSel;
out.m_q = m_q;
out.q = qWord;
out.q12 = q12;
out.nInstMin = nMin;
out.nInstMax = nMax;
out.nDacClip = nClip;
out.peEnd = pe;
out.vEnd = v;
out.mode = mode;
end
