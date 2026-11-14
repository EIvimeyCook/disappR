# Published figures against disappR (visual check)

Each sheet puts a published panel (left) beside the same analysis in the app (right): the bundled example, its
GAMM settings, step 7's trajectory plot (data shown as raw records, or no data where the paper shows none) and, for
Tsui et al., the Ageing rate plot. Built from app screenshots and the papers' figures (Tsui et al. 2026, Evol. Lett.,
CC BY 4.0, figures from the Europe PMC package; Sanghvi et al. 2021, Evol. Ecol., figures from the article PDF).

| Sheet | Panels |
|---|---|
| `tsui_fig1_vs_app_1..3.png` | Tsui Figure 1 (trajectories), 12 of 13 traits (telomere data not deposited) |
| `tsui_fig2_vs_app_1..3.png` | Tsui Figure 2 (ageing rates, first derivatives) |
| `sanghvi2021_vs_app.png` | Sanghvi Figures 2, 3A and 3C |

Expected differences:

- Intervals. The app draws simultaneous 95% intervals for the whole curve; Tsui et al. drew standard errors (Figure 1)
  and pointwise 95% intervals (Figure 2), which are narrower.
- Reference values. The app holds covariates and lifespan proxies at their means over individuals and averages factor
  covariates; papers often predict at other values. This shifts curves up or down without changing their shape: Tsui
  survival (app about 0.9, the observed survival of these birds; paper about 0.5) and reproductive success (app peak
  about 0.7, with last-reproduction age at its mean over birds; paper about 1.6).
- Families. Sanghvi Figure 3A used an ordinal (ocat) family; the app treats flight ability as continuous. Fecundity
  is Gaussian, as the published comparison appears to have been.
- Group labels. Tsui's sex is coded 0/1 in the deposit; the paper colours females red and males blue.
- Tsui survival leaves out the 76 records whose previous-year status is 'dead' (complete separation).

Shapes, peaks and the ordering of groups agree in every panel.
