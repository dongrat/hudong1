# AK2401A DMR two-point TX frequency-synthesis simulation

MATLAB/Octave simulation of **AK2401A-style two-point TX modulation** for **DMR 4FSK**: random symbols → RRC → Path A (ΣΔ PLL / frequency offset) + Path B (DAC → VCO) → combined frequency trajectory.

## How to run

```matlab
cd docs/ak2401a-dmr-twopoint-sim   % or absolute path to this folder
run_ak2401a_dmr_twopoint
```

- **MATLAB**: R2016b+ recommended (no toolboxes required).
- **Octave**: 8.x with `octave-signal` (`pkg load signal`) for `butter`; local fallbacks exist if missing.
- Headless: figures are also written under `figures/*.png` when `P.saveFigs = true`.

## Architecture

```
Random DMR 4FSK symbols  {-3,-1,+1,+3}
        │
        ▼
   RRC pulse shaping (α=0.2)  →  ideal Δf(t) [Hz]
        │
        ├──────────────► Path A command = LPF_loop(Δf)
        │                     │
        │                     ▼
        │              OFST quantize (18-bit)
        │                     │
        │                     ▼
        │              ΣΔ Frac-N + PLL closed-loop dynamics
        │                     │
        │                     ▼
        │                 f_PLL(t)
        │
        └──────────────► Path B command = Δf − LPF_loop(Δf)  (HF residual)
                              │
                              ▼
                       12-bit DAC + SMF (~20 kHz)
                              │
                              ▼
                       VCO: Δf = Km · v_SMF
                              │
                              ▼
                          f_VCO(t)

Combined:  f_comb(t) = f_PLL(t) + f_VCO(t)
RF model:  baseband FM envelope  exp(j 2π ∫ f_comb dt)  (carrier implied)
```

AK2401A maps to this as:

| Path | Chip block | Datasheet hook |
|------|------------|----------------|
| A | Frequency Offset **OFST1** → 18-bit ΔΣ Frac-N | Digital FM/FSK via offset; `f = fPFD·(INT+(FRAC+OFST)/MOD)` |
| B | **12-bit DAC + SMF → SMFOUT** | Connected to external VCO modulation / tune pin (“to VCO modulation”) |

The IC does not brand “two-point”; the sim uses the two TX paths the datasheet actually provides.

## Key parameters

| Parameter | Value | Origin |
|-----------|-------|--------|
| Modulation | DMR **4FSK** | ETSI TS 102 361; datasheet lists 4FSK |
| Symbol rate | **4800 baud** | DMR |
| Deviations | ±648 / **±1944 Hz** | ETSI (±1 / ±3) |
| RRC | **α = 0.2**, span 8 | DMR practice |
| TCXO / PFD | 19.2 MHz / **2.4 MHz** | Datasheet recommend / examples |
| MOD | **262143** | Datasheet “normally max” |
| Carrier | **450 MHz** | Assumed UHF DMR (external VCO) |
| PLL loop BW | **~800 Hz** | Assumed (external loop filter) |
| Kv (coarse tune) | **20 MHz/V** | Assumed external VCO (PLL) |
| Km (mod port) | **8 kHz/V** | Assumed SMFOUT→VCO mod sensitivity |
| DAC | **12-bit @ 96 kHz** | Datasheet max 200 kHz; S/N test uses 96 kHz |
| SMF | **fc ≈ 20 kHz** | Datasheet (−4 dB @ 20 kHz) |
| ICP | ≈ **640 µA** @ RBIAS2=27 kΩ, CPFINE=7 | `2160/R × (CPFINE+1)` |

Full extraction: `../../internal/ak2401a-sim/datasheet-params.md`.

## Files

| File | Role |
|------|------|
| `run_ak2401a_dmr_twopoint.m` | **Main entry** |
| `params_ak2401a_dmr.m` | Parameter struct |
| `gen_dmr_rrc_baseband.m` | Random symbols + RRC → Hz |
| `rrc_filter.m` | RRC FIR |
| `twopoint_split.m` | Complementary A/B commands |
| `model_ak2401a_sd_pll.m` | OFST + ΣΔ + PLL LPF |
| `model_ak2401a_vco_dac.m` | DAC + SMF + Kv |
| `plot_sim_results.m` | Plots |
| `figures/` | PNG outputs |
| `last_run_results.mat` | Last run data |

## What each plot shows

1. **01_symbols_rrc** — Random 4FSK symbols and RRC-shaped ideal deviation.
2. **02_twopoint_paths** — Path A/B commands; PLL output; DAC/VCO output.
3. **03_combined_vs_ideal** — Combined Δf vs ideal; residual error.
4. **04_spectra** — PSD of ideal / combined / Path A / Path B deviation.
5. **05_fm_envelope_spectrum** — Equivalent FM complex-envelope spectrum.
6. **06_architecture_notes** — Parameter summary text figure.

## Assumptions (need user input only if you want a board-matched model)

1. **External VCO Kv**, **modulation-port Km**, and **loop-filter BW** are not in the AK2401A datasheet — set `P.Kv_HzPerV` / `P.Km_HzPerV` / `P.fLoop_Hz` to your board.
2. Two-point **gain trim** (`P.pathBGain`) is ideal (=1); real boards need calibration.
3. ΣΔ is a **baseband-rate** 2nd-order error-feedback sketch, not a full PFD-rate MASH bit-stream.
4. OFST is updated every baseband sample (96 kHz model); hardware SPI/register timing may be slower.
5. Full GHz RF carrier is **not** sampled; combined Δf + complex envelope is the equivalent model.
6. TX-only RRC has ISI, so instantaneous peaks can exceed ±1944 Hz; combined residual vs ideal is typically ~150 Hz RMS with default knobs (OFST LSB, ΣΔ sketch, SMF, tracking lag).
