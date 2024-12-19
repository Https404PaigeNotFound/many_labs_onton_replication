 % alpha (typically 8-12 hz)

%{

WE ARE NOW DOING A ANOVA IV:2(condition)x3(memory load) DV: Baseline-normalised alpha power from the alpha component. 

Process for H3: Alpha Power for Ignore vs Memorise
Identify Alpha Components:

Use ICA (already performed in preprocessing) to identify components with significant alpha power (8–12 Hz) during the task. This requires:
Examining power spectra of components to find those dominated by alpha power.
Confirming dipole locations for these components (via DIPFIT) near expected cortical regions (e.g., occipital or parietal regions for visual alpha, or frontal regions for task-related alpha).
These identified alpha components are the input for further analysis.
Compute ERSPs for Alpha Components:

For each identified alpha component, compute event-related spectral perturbations (ERSPs) for the Ignore and Memorise epochs separately.
Baseline normalisation: Subtract the mean log power during the fixation baseline from the log power during task epochs to normalise.
Statistical Analysis:

Compare the ERSP-derived alpha power between the Ignore and Memorise conditions.
Use either a paired t-test or a within-subject ANOVA depending on how you want to model the data:
t-test: If only one Ignore and one Memorise condition per subject (average alpha power across trials for each condition).
ANOVA: If you include additional factors, such as memory load (3, 5, or 7 letters) or time windows, making it a Condition (Ignore vs Memorise) × Load design.
Analyse differences in alpha power during Ignore vs Memorise epochs to test the hypothesis that alpha is involved in filtering irrelevant information.
What Does Onton et al. Do?
Onton et al. (2005) do not explicitly analyse alpha power in Ignore vs Memorise conditions.
They note the presence of alpha band activity but focus mainly on theta and beta.
Their exploratory findings suggest mixed effects for alpha power, with no formal hypothesis-driven analysis.
What Does Our Write-Up Say?
In your current write-up:

You plan to identify alpha components in a similar way to theta components but without detailing the specific selection criteria.
The analysis involves computing ERSPs for Ignore and Memorise conditions, with differences assessed statistically (likely paired t-tests or ANOVA).
Comments indicate concern about the lack of clarity in:
How alpha components are selected.
The justification for using ICA-space ERSPs for alpha, given it is not standard.
What Do the Comments and Feedback Say?
Alpha Component Selection:

Acknowledge that selecting alpha components in ICA-space is less standard than channel-based analysis.
Suggest clearer criteria for identifying alpha components (e.g., by power spectra or task-specific activation).
Suggest averaging across multiple alpha components if more than one is identified, with explicit criteria for how these are selected.
Statistical Analysis:

Feedback suggests ensuring the analysis reflects the hypothesis explicitly and uses appropriate methods.
For Ignore vs Memorise, the feedback indicates that a paired t-test may suffice unless you want to include additional factors like memory load or time.
Suggested Process for H3
1. Identify Alpha Components:
Use ICA output to identify components with:
Strong alpha power (8–12 Hz) in the power spectra.
Dipole locations in plausible cortical regions for alpha (frontal, parietal, occipital).
If multiple alpha components are identified:
Average across all alpha components for each subject.
Ensure averaging is justified and reported transparently.
2. Compute ERSPs:
Compute baseline-normalised ERSPs for alpha components during Ignore and Memorise epochs.
Use the same method as for theta (subtract fixation baseline).
3. Statistical Analysis:
Test for significant differences in alpha power between Ignore and Memorise:
Paired t-test: Compare average alpha power for Ignore vs Memorise conditions (across all trials per subject).
ANOVA: If memory load or time windows are included, use a Condition (Ignore vs Memorise) × Load (3, 5, 7 letters) within-subject design.
Conclusion
Process: Identify alpha components (via ICA), compute ERSPs for Ignore and Memorise epochs, and analyse using t-tests or ANOVA.
Statistical Test:
Paired t-test: Simpler, for single condition comparisons.
ANOVA: More flexible, allowing inclusion of additional factors like memory load.
What to Clarify:
Explicit criteria for alpha component selection.
Justify using ICA-space ERSPs (as this is non-standard for alpha analysis).
%}