# miFISH
Simultaneous visualization of DNA loci in single cells by combinatorial color iFISH

miFISH is an extension of the recently developed iFISH method (Gelali et al, Nature Communications 2019), which allows to localize multiple chromosome loci in situ with high accuracy through combinatorial labeling. The datasets are valuable for researchers interested in nuclear architecture on the micrometer scale and we make available the codes necessary to analyse the datasets and to decode the corresponding loci for each dot.

The dots coordinates for each FISH dot is extracted with DOTTER as with iFISH method. To decode the signals, we present the following pipelines that were described in detail in Mota et al, Scientific Data in revision. We recomment the miFISH pipeline with all filters included because it provides the best correlation with HiC data. 


<br>

### miFISH pipeline with all filters (best correlation with HiC)
Find the script for this analysis in **_miFISH experiment > Signal decoding > with Filters_**.
Extract the 2 scripts to a working folder where the **_datasets folder_** should be in the same directory as the working script.
The script **_main_miFISH_** is ready to be run by the user. To specify the dataset replicate to be run, change the number in line 6 from "1" to "2".

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
Extract the 2 scripts to a working folder where the **_datasets folder_** should be in the same directory as the working script.
The script **_main_miFISH_** is ready to be run by the user. To specify the dataset replicate to be run, change the number in line 6 from "1" to "2".

- detect 10 signals for each fluorescent dye (except for single-color probes which are immediately assigned the brighest of the channel);
- compute all pairwise distances in between all the signals for all the channels, to match signals coming from the same dual-color probe;
- examine dots with pairwise distances lower than the defined cut-off distance;
The same dots can be allocated to other probes, which is in contrary to the miFISH pipeline with filters.

<br>

### miFISH stochastic allocation of dots
Find the script for this analysis in **_miFISH experiment > Signal decoding > Stochastic_**.
Extract the 2 scripts to a working folder where the **_datasets folder_** should be in the same directory as the working script.
The script **_main_miFISH_** is ready to be run by the user. To specify the dataset replicate to be run, change the number in line 6 from "1" to "2".

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

## miFISH vs HiC
Find the script for this analysis in **_miFISH experiment > miFISH vs HiC_** folder.
Extract the script to a working folder where the **_datasets folder_** should be in the same directory as the working script.
Extract and run one of the previous miFISH pipeline with/without all filters or stochastic allocation. The "Coordinates" and "Lamina" variables are re-used by the HiC script.
Lastly, run the script **_miFISH vs HiC_**.

<br>

## Pairwise distances, Lamina distances and Compartments
Extract the script for this analysis in **_miFISH experiment > Lamina distances and Compartments_** folder.
Additionally, extract the script in **_miFISH experiment > Signal decoding > with Filters_**.
In line 165 of **_main_miFISH_** remove the "%" in "% laminaDISTeigen(Coordinates, Lamina)".
Run the analysis with all 3 scripts in the same working directory with the folder **_datasets_**.


Figures:
- All pairwise distances combination of probes with the corresponding eigenvector combination;
- Boxplot representation for each eigenvector combination (A-A) (A-B) (B-B);
- Taking all pairwise distances for the 120 probes combinations, measure the median for each 120 combination and make a boxplot that divides the 120 distances for each eigenvector combination (A-A) (A-B) (B-B) (or no compartment asigned);
- All pairwise distances combination of probes (120 in total) represented according to the lamina distance differences in between the probes;
- Representation of the median of lamina distance normalised for each probe and the color code displays the compartment classification;
- Boxplot for all lamina distances measure for each probe;
- Representation of the Eigenvector distribution for all chr2, the probes are positioned at their corresponding genomic coordinate and the number indexed is the median lamina distance;
- All 3 and 20 distances grouped and show the behaviour of each 3Mbp and 20 Mbp.

<br>

## Full Width at Half Maximum
Find the script for this analysis in **_miFISH experiment > FWHM_** folder.
Extract the script to a working folder where the **_datasets folder_** should be in the same directory as the working script.
The script **_main_no_filter_fwhm_** is ready to be run by the user. To specify the dataset replicate to be run, change the number in line 6 from "1" to "2".

<br>

## Structures curvature
Find the script for this analysis in **_miFISH experiment > curvature_** folder.
Extract the 2 scripts to a working folder where the **_datasets folder_** should be in the same directory as the working script.
The script **_main_curvature_** is ready to be run by the user. To specify the dataset replicate to be run, change the number in line 6 from "1" to "2".

<br>

# Extra experiments of miFISH with one dual-color
## Cut-off for different channels matching
The co-localization calling for dual-color probes might be challenging because the pairwise distance in between the two color dots is not 0 due to residual chromatic aberration that was not completely corrected and to other technical factors. Therefore another set of experiments were performed to measure this technical error and identify a cut-off for co-localization calling. 

We performed two sets of experiments, in which we used only one dual-color probe alongside two singly labelled probes targeting different loci on chr2. We tested the shift for two different dual-color probes, where one of them can be imaged using the same multi-band dichroic mirror, while the other requires a change of the cube holding dichroic mirrors. The median distance between two signals coming from the same probe labelled with ATTO 542 (tmr) and ATTO 647N (a647)—for which dichroic mirrors are placed in the same cube—was lower than 0.25 µm. In contrast, the median distance between two signals coming from the same probe labelled with Alexa Fluor 488 (a488) and Alexa Fluor 594 (a549)—for which dichroic mirrors are placed in different cubes and hence more mechanical movement is needed to image both—increased to 0.55. Therefore, in further analyses, we retained only signals with pairwise distance lower than 0.55 µm for experiments in which dichroic mirrors were placed in different cubes, while for experiments in which dichroic mirrors shared the same cube, we set the threshold to 0.25 µm.

Find the script for this analysis in **_miFISH colour overlap distance experiment_** folder.
Extract the script to a working folder where the **_datasets folder_** should be in the same directory as the working script.
The script **_main_single_probes_** is ready to be run by the user. To specify one of the overlap color tests to be run, change the dyes combination written in line 5 from "AT647N-AT542" to "AF594-AF488".

<br>

# Extra experiments of iFISH for chr1, chr2 and chr10
We performed additional iFISH experiments using a set of 60 singly labelled iFISH probes targeting multiple evenly spaced loci along chr1, 2 and 10, which we retrieved from our iFISH probe repository (Gelali et al Nat. Comm. 2019). We then calculated all the pairwise distances either between probes labelled with the same color or between probes labeled with two different colors and examined how the pairwise distances change with increasing genomic distances between the probes.

Find the script for this analysis in **_iFISH chr1 chr2 chr10_** folder.
Extract the script to a working folder where the **_datasets folder_** should be in the same directory as the working script.
The script **_main_chr1_2_10_single_colour_** is ready to be run by the user. To specify the chromosome, change the number in line 6 from "2" to "1" or "10". To specify the dataset replicate to be run, change the number in line 10 from "_rep1" to "_rep2" or for chr1 is empty "" (no replicate).
