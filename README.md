# miFISH
Simultaneous visualization of DNA loci in single cells by combinatorial color iFISH

miFISH is an extension of the recently developed iFISH method (Gelali et al, Nature Communications 2019), which allows to localize multiple chromosome loci in situ with high accuracy through combinatorial labeling. The datasets are valuable for researchers interested in nuclear architecture on the micrometer scale and we make available the codes necessary to analyse the datasets and to decode the corresponding loci for each dot.

The dots coordinates for each FISH dot is extracted with DOTTER as with iFISH method. To decode the signals, we present the following pipelines that were described in detail in Mota et al, Scientific Data in revision. We recomment the miFISH pipeline with all filters included because it provides the best correlation with HiC data. 

<br>

### miFISH pipeline with all filters (best correlation with HiC)
Find the script for this analysis in **_miFISH experiment > Signal decoding > with Filters_**.

- detect 10 signals for each fluorescent dye (except for Alexa Fluor 790 (ir800) channel is expected 2 dots);
- compute all pairwise distances in between all the signals for all the channels, to match signals coming from the same dual-color probe;
- examine dots with pairwise distances lower than the defined cut-off distance;
- prioritize dual-color probes that only had one possibility to occur;
- implemente a cut-off value to filter out nearby dots of the same color, given that the signal originating from a single FISH probe might split into two or more dots;
- filter out dots positioned more than 5 µm away from all the other dots in the same cluster;
- assign single-color probes to the brightest dots of the cluster that had not yet been selected. 
 
<br>

### miFISH pipeline simplified parameters
Find the script for this analysis in **_miFISH experiment > Signal decoding > without Filters_**.

- detect 10 signals for each fluorescent dye (except for single-color probes which are immediately assigned the brighest of the channel);
- compute all pairwise distances in between all the signals for all the channels, to match signals coming from the same dual-color probe;
- examine dots with pairwise distances lower than the defined cut-off distance;
The same dots can be allocated to other probes, which is in contrary to the miFISH pipeline with filters.

<br>

### miFISH stochastic allocation of dots
Find the script for this analysis in **_miFISH experiment > Signal decoding > Stochastic_**.

- assign an aleatory dot to a single-color probe or a dual-color probe that corresponds to one of the required channels;
- in case of a dual-color find the other color dot that is positioned closer to the aleatory dot picked in first place; 
- in case the pairwise distance is higher than 0.55 µm, do not assign any dot and start again with another aleatory dot;
- measure for every iteration out of 100, the best combination of dots by applying an objective function:

  1 - Brightness intensity, *e_brightness* <br>
  2 - number of repeated dots in the set, *e_repeat_dots* <br>
  3 - overlap distance (0-0.55 µm), *e_probe_intra_dist* <br>
  4 - number of missing probes in the set, *e_dot* <br>
  5 - closeness to other dots from same channel (penalty to distances < 0.25µm), *e_merge_dot*

  Total error = e_brightness + 3 x e_repeat_dots + 10 x e_probe_intra_dist + 10 x e_dot + 10 x e_merge_dot

- for each allele, the procedure is re-started 20 times and only saves the set with the best overall score of all probes defined by the objective function.<br>

<br>

## Pairwise distances
After choosing one of the previous decoding methods and running the script, extract the script for this analysis in **_miFISH experiment > Pairwise distances_** folder.

<br>

## miFISH vs HiC
After choosing one of the previous decoding methods and running the script, extract the script for this analysis in **_miFISH experiment > miFISH vs HiC_** folder.

<br>

## Lamina distances and Compartments
After choosing one of the previous decoding methods and running the script, extract the script for this analysis in **_miFISH experiment > Lamina distances and Compartments_** folder.

Figures:
- All pairwise distances combination of probes with the corresponding eigenvector combination;
- Boxplot representation for each eigenvector combination (A-A) (A-B) (B-B);
- Taking all pairwise distances for the 120 probes combinations, measure the median for each 120 combination and make a boxplot that divides the 120 distances for each eigenvector combination (A-A) (A-B) (B-B) (or no compartment asigned);
- All pairwise distances combination of probes (120 in total) represented according to the lamina distance differences in between the probes;
- Representation of the median of lamina distance normalised for each probe and the color code displays the compartment classification;
- Boxplot for all lamina distances measure for each probe;
- Representation of the Eigenvector distribution for all chr2, the probes are positioned at their corresponding genomic coordinate and the number indexed is the median lamina distance.

<br>

## Full Weight Height Maximum
Find the script for this analysis in **_miFISH experiment > FWHM_** folder.

<br>

# Extra experiments of miFISH with one dual-color
## Cut-off for different channels matching
The co-localization calling for dual-color probes might be challenging because the pairwise distance in between the two color dots is not 0 due to residual chromatic aberration that was not completely corrected and to other technical factors. Therefore another set of experiments were performed to measure this technical error and identify a cut-off for co-localization calling. 

We performed two sets of experiments, in which we used only one dual-color probe alongside two singly labelled probes targeting different loci on chr2. We tested the shift for two different dual-color probes, where one of them can be imaged using the same multi-band dichroic mirror, while the other requires a change of the cube holding dichroic mirrors. The median distance between two signals coming from the same probe labelled with ATTO 542 (tmr) and ATTO 647N (a647)—for which dichroic mirrors are placed in the same cube—was lower than 0.25 µm. In contrast, the median distance between two signals coming from the same probe labelled with Alexa Fluor 488 (a488) and Alexa Fluor 594 (a549)—for which dichroic mirrors are placed in different cubes and hence more mechanical movement is needed to image both—increased to 0.55. Therefore, in further analyses, we retained only signals with pairwise distance lower than 0.55 µm for experiments in which dichroic mirrors were placed in different cubes, while for experiments in which dichroic mirrors shared the same cube, we set the threshold to 0.25 µm.

Find the script for this analysis in **_miFISH colour overlap distance experiment_** folder.

<br>

# Extra experiments of iFISH for chr1, chr2 and chr10
We performed additional iFISH experiments using a set of 60 singly labelled iFISH probes targeting multiple evenly spaced loci along chr1, 2 and 10, which we retrieved from our iFISH probe repository (Gelali et al Nat. Comm. 2019). We then calculated all the pairwise distances either between probes labelled with the same color or between probes labeled with two different colors and examined how the pairwise distances change with increasing genomic distances between the probes.

Find the script for this analysis in **_iFISH chr1 chr2 chr10_** folder.
