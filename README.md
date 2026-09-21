# Spatial Transport Identifiability in Desmoplastic Tumours: When a Lumped Burden ODE Cannot Represent a Fibrotic Delivery Barrier

For stroma-rich carcinomas, which observables (interstitial pressure, collagen anisotropy, perfusion maps) are required before a zero-dimensional burden ODE is even structurally capable of representing delivery failure, versus remaining an unidentified lumped sink?

**Thesis #11** (series label NP-02). Computational research, set out in Nile University B.Sc. chapter order for handoff.

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**GitHub:** https://github.com/cloudynirvana  
**Date:** 21 September 2026

The calculations use one three-shell radial compartment reduction of a filtration-diffusion balance, and one lumped burden ODE with a single delivery coefficient. Barrier coordinates are a hydraulic factor, a collagen fraction, a circumferential anisotropy, and an interstitial pressure. They enter the shells only through two conductances. Schedule AUC records the time-integral of the spatial mean. Schedule L records that mean at 24 times. Schedule P records the three shells. Schedule PI adds one pressure reading. Schedule PIA adds one anisotropy reading.

On this generator the integral has Fisher rank 1 of 4. The mean time course has structural rank 2 and practical rank 1: the second principal relative standard error is 3.21. A noiseless fit of the lumped series, started away from the truth, stops on a different barrier vector whose shell trajectories match the truth to solver tolerance. Perfusion raises practical rank to 2 of 4, which is the identifiable subset carried by the two conductances. Pressure with perfusion raises rank to 3. Pressure, perfusion, and anisotropy raise rank to 4, and a noiseless fit from the same distant start returns the truth.

A second barrier, retuned so that the mean stays close, drops the core concentration at 8 h from 0.00947 to 0.000816. The lumped ODE, having one state, cannot represent that split.

Draws are a synthetic surrogate. Chapter Four is not a perfusion scan and not a collagen measurement.

This is research only. It is not a medical device, not clinical decision support, not a dose, and not a cure. No document DOI is registered.

See [DISCLAIMER.md](DISCLAIMER.md). The manuscript is [THESIS.md](THESIS.md).

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Manuscript (Chapters 1 to 5, Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata, no document DOI |
| `DISCLAIMER.md` | Research-only boundary |
| `sim/transport_identifiability.py` | Seeded compartment, lumped, and Fisher sketches (seed 20260921) |
| `sim/results.json` | Numbers cited in Chapter Four |
| `sim/figures/` | Trajectories, core comparison, spectra, and level sets |

## Reproduce

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/transport_identifiability.py
python3 build_pdf.py
```

NumPy, SciPy and Matplotlib are required for the sketches. The PDF step also needs the `markdown` and `weasyprint` packages. Regenerating the script rewrites `sim/results.json` and `sim/figures/`.

## Cite

Ogbonna KE. Spatial transport identifiability in desmoplastic tumours: when a lumped burden ODE cannot represent a fibrotic delivery barrier [Internet]. Thesis #11 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-11-desmoplastic-transport-identifiability

Machine-readable fields are in `CITATION.cff`. Add a document DOI there only after one exists.

Hub index, for cataloguing only: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

## Licence

Text and sketch code are MIT, with attribution. Computational research only.
