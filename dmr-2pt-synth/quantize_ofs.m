function [mq, ofs, q12] = quantize_ofs(P, m)
%QUANTIZE_OFS  Shared 12-bit quantizer, applied before either TX path.
%  Full-scale span matches the signed 16-bit OFS word. Only the upper 12
%  bits are kept, so the RF step is 16 times the OFS LSB.
%    q16 = round(m * OFS_DEN / fPFD)
%    q12 = clamp(round(q16 / 16), -2048, 2047)
%    m_q = q12 * 16 * fPFD / OFS_DEN
%    Lch OFS_MDLT = q12 * 16
%  Hch is computed only from m_q.

q16 = round(m .* (P.ofsDen / P.fPfd_Hz));
q12 = round(q16 / 16);
q12 = min(2047, max(-2048, q12));
ofs = q12 * 16;
mq = ofs .* (P.fPfd_Hz / P.ofsDen);
end
