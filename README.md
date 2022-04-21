# Azkanaz et al, 2022. Nature.

This repository contains the code and data used in Azkanaz et al. (Nature 2022). Simulation and fits were made using Python (custom pipelines and scipy.optimize library) and R (4.1.1) was used for RNAseq analysis and cell migration.

## Abstract
Morphology and functionality of the epithelial lining differ along the intestinal tract, yet tissue renewal at all sites is driven by stem cells at the base of crypts (Lopez-Garcia et al., 2010; Barker et al., 2007). Whether stem cell numbers and behaviour vary at different sites is unknown. Here, we show by intravital microscopy that despite similarities in  the number and distribution of proliferative cells with an Lgr5 signature, small intestinal (SI) crypts contain twice as many effective stem cells as large intestinal (LI) crypts. We find that, although passively displaced by a conveyor belt-like upward movement, SI cells positioned away from the crypt base can function as long-term effective stem cells due to Wnt-dependent retrograde cellular movement. By contrast, the near absence of retrograde movement in the LI restricts cell repositioning, leading to a reduced effective stem cell number. Moreover, upon suppression of the retrograde movement in the SI, the number of effective stem cells is reduced, and the rate of monoclonal conversion of crypts is accelerated. Together, these results show that effective stem cell number is determined by active retrograde movement, revealing a new channel of stem cell regulation that can be experimentally and pharmacologically manipulated.

## - Simulation analysis
Python scripts/ integration of mathematical solutions and fitting/ Simulation of recovery dynamics.

- Clone_survival_mathematical.py: Fits the kr/kd paramters from real data at the last time point using scipy.optimize.curve_fit. Computes the probability of clone retention in time according to the mathematical prediction provided in equation 1.11 of the supplementary note. It outputs the plot corresponding to Ext. Data Fig 9e.

- Comparison_Numerics_Theory.py: Compares the results obtained through the simulation of the SCB dynamics in 1D, provided in the file Modelling_data/Model_Predictions_General.txt, to the mathematical predictions provided by equation 1.11 of the supplementary note. It generates as output the plot corresponding to Ext. Data Fig. 9i.

- Deepmost_Cell_Ablation.py: Simulation of crypt recovery. Steps are described in the supplementary theory section 2.5. Outputs the plots corresponding to figure 4h,i.

- Migration_to_kr_header.py: Header containing the functions to run the above described scripts.

## - 2D simulation of crypts
Parameters of the simulations are found in "variables.cpp", main script is "main.cpp" and raw results outputted in "write_data.cpp". This compiles on standard linux by writing in the terminal:

```bash
./clean.sh; g++ -o sim main.cpp -lm  -lgsl -lgslcblas; ./sim
```

Raw output (statistics for 1000 independent repetitions of n clonal inductions, used to build up averages and confidence intervals) is stored in "surv.txt". Run afterwards the Python script:

```bash
python averaging.py
```

To give rise to average survival probabilities as a function of starting position and time  (day 2, day 3, day 4, dat 56). This can then be plotted using any software, for instance, by GNUplot, version 5.4 with the script:

```gnuplot
set style fill transparent solid 0.1 noborder; a=1
set xrange [0:3.5]
set yrange [0:1.1]
set xlabel "Position"
set ylabel "Survival probability"
set tics nomirror

plot "model_pos.txt" u ($0/a):12:17 with filledcurves lt 2 title "Day 2", "" u ($0/a):13:18 with filledcurves lt 3 title "Day 3", "" u ($0/a):14:19 with filledcurves lt 4 title "Day 4", "" u ($0/a):15:20 with filledcurves lt 7 title "Day 56"

replot "model_pos.txt" u ($0/a):2 w l lw 2 lt 2 notitle, "" u ($0/a):3 w l lw 2 lt 3 notitle, "" u ($0/a):4 w l lw 2 lt 4 notitle, "" u ($0/a):5 w l lw 2 lt 7 notitle
```
Simulation analysis and 2D simulation of crypts correspond to: figures 2d, 3g, 3h, 4c, 4d, 4h and extended figures 9e, 9f, 9g, 9h, 9i, 9j, 9k, 9l, 9m, 10c. 

## - RNAseq analysis
All code for RNAseq analysis was developed using DESeq2 in R version 4.1.1. All RNAseq data is publicly available through Gene Expression Omnibus (GEO): [GSE194250](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194250). The corresponding figures are: figure 1h, Extended figure 1c, 1d.

## - In vitro migration assay
The corresponding figures are: figures 3a, 3b and extended figure 7c

## - Cell tracking analysis
All code for cell tracking was developed on Python and Jupyter notebooks. All data is included. The corresponding figure is: figure 3d.

## - Data availability
All data is available, see Source data of the original manuscript. 

## Contact and support
- Maria Azkanaz. Netherlands Cancer Institute. m.azkanaz@nki.nl
- Jacco van Rheenen (Principal investigator). Netherlands Cancer Institute. j.v.rheenen@nki.nl

## References
Lopez-Garcia, C., Klein, A. M., Simons, B. D. & Winton, D. J. Intestinal stem cell replacement follows a pattern of neutral drift. Science 330, 822-825, doi:10.1126/science.1196236 (2010).

Barker, N. et al. Identification of stem cells in small intestine and colon by marker gene Lgr5. Nature 449, 1003-1007, doi:10.1038/nature06196 (2007).
set style fill transparent solid 0.1 noborder; a=1; set xrange [0:3.5]; set yrange [0:1.1]; set xlabel "Position"; set ylabel "Survival probability" ; set tics nomirror;
plot "model_pos.txt" u ($0/a):12:17 with filledcurves lt 2 title "Day 2", "" u ($0/a):13:18 with filledcurves lt 3 title "Day 3", "" u ($0/a):14:19 with filledcurves lt 4 title "Day 4", "" u ($0/a):15:20 with filledcurves lt 7 title "Day 56";
replot "model_pos.txt" u ($0/a):2 w l lw 2 lt 2 notitle, "" u ($0/a):3 w l lw 2 lt 3 notitle, "" u ($0/a):4 w l lw 2 lt 4 notitle, "" u ($0/a):5 w l lw 2 lt 7 notitle;