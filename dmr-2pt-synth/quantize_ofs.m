function [mq, q] = quantize_ofs(P, m)
%QUANTIZE_OFS  Signed 16-bit OFS quantizer, applied before either TX path.
%  LSB is the OFS frequency step at RF: fPFD / OFS_DEN.
%    q   = clamp(round(m * OFS_DEN / fPFD), -32768, 32767)
%    m_q = q * fPFD / OFS_DEN
%  Lch writes q straight into OFS_MDLT. Hch is computed only from m_q.

q = round(m .* (P.ofsDen / P.fPfd_Hz));
q = min(32767, max(-32768, q));
mq = q .* (P.fPfd_Hz / P.ofsDen);
end
