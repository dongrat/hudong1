function L = design_loop_filter(P)
%DESIGN_LOOP_FILTER  3rd-order passive-equivalent loop filter.
%  Crossover and phase margin follow Banerjee's type-II design, then a
%  high-frequency pole is added. Component values realize the datasheet's
%  ~1 kHz evaluation-board bandwidth; R and C themselves are not published.

wc = 2 * pi * P.fLoop_Hz;
pm = P.phaseMargin_deg * pi / 180;
% Lead network: geometric-mean placement gives the requested phase boost.
T1 = (1 / cos(pm) - tan(pm)) / wc;
T2 = 1 / (wc^2 * T1);

Kpd = P.Icp_A / (2 * pi);          % A/rad
Kv_rad = 2 * pi * P.Kvco_HzPerV;   % (rad/s)/V
N = P.fVco_Hz / P.fPfd_Hz;

C1 = (T1 / T2) * (Kpd * Kv_rad) / (wc^2 * N) ...
    * sqrt((1 + (wc * T2)^2) / (1 + (wc * T1)^2));
C2 = C1 * (T2 / T1 - 1);
R2 = T2 / C2;
Ctot = C1 + C2;

T3 = 1 / (2 * pi * P.fPole3_Hz);
% Equivalent series-R / shunt-C for the extra pole (not a board extract).
C3 = 1e-9;
R3 = T3 / C3;

L.T1 = T1;
L.T2 = T2;
L.T3 = T3;
L.C1_F = C1;
L.C2_F = C2;
L.Ctot_F = Ctot;
L.R2_ohm = R2;
L.C3_F = C3;
L.R3_ohm = R3;
L.N = N;
L.Kpd = Kpd;
L.Icp_A = P.Icp_A;
end
