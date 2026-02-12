# Experimental Design Document

## Study: Discrete Choice Experiment on Breast Cancer Screening Preferences

### Protocol: DCE-BCS-2026-001
### Version: 1.0 | Date: February 2026

---

## 1. Design Overview

### 1.1 Attributes and Levels

| Attribute | Levels | Coding |
|-----------|--------|--------|
| Screening method | Mammography, MRI, Ultrasound | Effects-coded (2 dummy variables) |
| Screening frequency | Every year, Every 2 years, Every 3 years | Effects-coded (2 dummy variables) |
| Out-of-pocket cost | $0, $50, $150, $300 | Continuous (linear) |
| Sensitivity (detection rate) | 70%, 85%, 95% | Continuous (linear) |
| False positive rate | 5%, 10%, 15% | Continuous (linear) |
| Waiting time for results | 1 day, 1 week, 3 weeks | Effects-coded (2 dummy variables) |
| Pain/discomfort | None, Mild, Moderate | Effects-coded (2 dummy variables) |

**Total attribute levels:** 3 + 3 + 4 + 3 + 3 + 3 + 3 = 22 levels across 7 attributes

### 1.2 Full Factorial Size

The full factorial design would produce:

3 x 3 x 4 x 3 x 3 x 3 x 3 = **2,916 possible profiles**

For a paired comparison (2 alternatives per choice set), the number of possible choice sets is:

C(2916, 2) = 4,252,170 choice sets

This is far too large to administer; hence, a fractional factorial / optimal design approach is required.

---

## 2. Design Strategy

### 2.1 Approach: D-optimal Fractional Factorial Design

We adopt a **D-optimal design** approach rather than an orthogonal main-effects plan (OMEP) for the following reasons:

1. **Flexibility**: D-optimal designs can accommodate any number of attributes and levels, including unbalanced designs
2. **Efficiency**: They maximize the determinant of the Fisher information matrix, producing more precise parameter estimates
3. **Practical constraints**: They allow incorporation of constraints (e.g., avoiding implausible combinations)
4. **Statistical power**: D-optimal designs outperform random or orthogonal designs in terms of estimation efficiency for the same number of choice sets

### 2.2 Alternative Considered: Orthogonal Main Effects Plan

An orthogonal fractional factorial design was considered but rejected because:
- Orthogonality is defined at the profile level, not the choice-set level
- It does not account for the pairing of profiles into choice sets
- It cannot guarantee level balance across alternatives within choice sets
- It produces lower D-efficiency compared to algorithmic optimal designs

### 2.3 Design Properties

The design targets the following properties:

| Property | Target |
|----------|--------|
| **D-efficiency** | >= 90% relative to theoretical optimum |
| **A-efficiency** | >= 85% |
| **Level balance** | Each level appears approximately equally often across all choice sets |
| **Orthogonality** | Minimal correlation between attributes within alternatives |
| **Minimal overlap** | Attribute levels differ between alternatives within a choice set to maximize information |
| **Utility balance** | Alternatives are approximately equally attractive (based on prior estimates) |

---

## 3. Design Generation

### 3.1 Software Recommendation

**Primary: Ngene (version 1.3+)**

Ngene is the recommended software for generating the experimental design because:
- Purpose-built for stated choice experiments
- Supports Bayesian D-optimal designs with informative priors
- Can handle MNL, mixed logit, and latent class model specifications
- Provides comprehensive design diagnostics (D-error, S-estimate, efficiency measures)
- Widely used in health economics DCE literature

**Alternative options:**
- **R (ideXact or AlgDesign packages)**: Open-source alternative; `ideact` supports D-optimal designs for DCEs; `AlgDesign` for general optimal design
- **JMP (Choice Design platform)**: Commercial; good GUI but less flexible than Ngene for Bayesian designs
- **Sawtooth Software (Lighthouse Studio)**: Integrated survey + design platform; good for implementation but less control over design properties

### 3.2 Ngene Design Syntax (Draft)

```
Design
;alts = alt1, alt2
;rows = 48
;block = 4
;eff = (mnl,d)
;model:
U(alt1) = b0
  + b1.effects[0.1,-0.1] * method[1,2,3]
  + b2.effects[0.05,-0.05] * frequency[1,2,3]
  + b3[-0.005] * cost[0,50,150,300]
  + b4[0.03] * sensitivity[70,85,95]
  + b5[-0.02] * falsepositive[5,10,15]
  + b6.effects[-0.1,0.05] * waittime[1,7,21]
  + b7.effects[0.15,-0.05] * pain[0,1,2] /
U(alt2) = b0
  + b1 * method
  + b2 * frequency
  + b3 * cost
  + b4 * sensitivity
  + b5 * falsepositive
  + b6 * waittime
  + b7 * pain
$
```

### 3.3 Prior Parameter Values

Prior parameter estimates for Bayesian D-optimal design are drawn from:
1. Published DCE studies on cancer screening preferences
2. Expert elicitation from clinical collaborators
3. Results from pilot study (to be conducted)

| Parameter | Prior mean | Prior SD | Source |
|-----------|-----------|----------|--------|
| Method (MRI vs. ref) | 0.10 | 0.15 | Literature |
| Method (Ultrasound vs. ref) | -0.10 | 0.15 | Literature |
| Frequency (2yr vs. ref) | 0.05 | 0.10 | Expert |
| Frequency (3yr vs. ref) | -0.05 | 0.10 | Expert |
| Cost (per $1) | -0.005 | 0.003 | Literature |
| Sensitivity (per %) | 0.03 | 0.02 | Literature |
| False positive (per %) | -0.02 | 0.015 | Literature |
| Wait time (1 week vs. ref) | -0.10 | 0.10 | Expert |
| Wait time (3 weeks vs. ref) | 0.05 | 0.10 | Expert |
| Pain (Mild vs. ref) | 0.15 | 0.10 | Literature |
| Pain (Moderate vs. ref) | -0.05 | 0.10 | Literature |

---

## 4. Choice Set Configuration

### 4.1 Number of Choice Sets

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| **Total choice sets** | 48 | Sufficient to estimate all main effects and selected interactions |
| **Blocks** | 4 | Each respondent sees 12 choice sets (manageable cognitive burden) |
| **Choice sets per respondent** | 12 | Within the recommended range of 8-16 for health DCEs |
| **Alternatives per choice set** | 2 + opt-out | Forced choice with status quo / opt-out option |
| **Survey versions** | 4 | One per block; respondents randomly assigned |

### 4.2 Blocking Strategy

- 48 total choice sets divided into 4 blocks of 12
- Blocks are orthogonal to attribute levels (blocking variable is uncorrelated with design columns)
- Respondents are randomly assigned to one block
- Each block is a self-contained, efficient sub-design

### 4.3 Opt-out Option

An opt-out / status quo option ("I would choose neither program") is included because:
1. It reflects real-world decision making (women can choose not to screen)
2. It avoids forced-choice bias
3. It allows estimation of participation rates and welfare measures
4. Required for valid willingness-to-pay calculations

The opt-out is coded as a third alternative with all attributes set to zero (no screening).

---

## 5. Design Quality Checks

### 5.1 Dominance Test

Choice set 8 is designed as a **dominance test** where Program A is strictly superior to Program B on all attributes. Respondents who fail to choose the dominant option may be:
- Not paying attention
- Responding randomly
- Failing to understand the task

**Decision rule:** Flag respondents who fail the dominance test. Conduct sensitivity analysis with and without these observations.

### 5.2 Stability Test (Test-Retest)

One choice set is repeated later in the survey (choice set 12 mirrors choice set 3 with alternatives swapped). This provides a measure of response consistency.

**Decision rule:** Respondents with inconsistent responses are flagged but not automatically excluded. Conduct sensitivity analysis.

### 5.3 Design Diagnostics to Report

- D-error (lower is better)
- A-error
- B-estimate (sample size needed per Bayesian approach)
- S-estimate (sample size estimate from design)
- Attribute level balance statistics
- Pairwise correlation matrix of attribute levels
- Condition number of the design matrix

---

## 6. Design Efficiency Metrics

### 6.1 D-efficiency

D-efficiency is defined as:

$$D_{eff} = \left(\frac{|\Omega^{-1}|^{1/K}}{N}\right)^{-1}$$

Where:
- $\Omega$ is the variance-covariance matrix of parameter estimates
- $K$ is the number of parameters
- $N$ is the number of choice sets

**Target:** D-efficiency >= 90%

### 6.2 Relative D-efficiency

We compute relative D-efficiency compared to a theoretical benchmark (full factorial enumeration):

$$D_{rel} = \frac{D_{eff}(\text{fractional})}{D_{eff}(\text{full factorial})} \times 100\%$$

### 6.3 S-estimate

The S-estimate indicates the minimum sample size required to obtain statistically significant parameter estimates at the 95% confidence level:

$$S = \frac{K}{D_{eff} \times T \times J}$$

Where T = number of choice tasks per respondent and J = number of alternatives.

---

## 7. Model Specification

### 7.1 Primary Model: Mixed Logit (Random Parameters)

The design supports estimation of a mixed logit model:

$$U_{njt} = \beta_n' x_{njt} + \varepsilon_{njt}$$

Where:
- $\beta_n$ ~ N($\mu$, $\Sigma$) captures preference heterogeneity
- $x_{njt}$ is the attribute vector for alternative $j$ in choice set $t$ for individual $n$
- $\varepsilon_{njt}$ is iid Type I Extreme Value

### 7.2 Secondary Models

1. **Conditional logit (MNL)**: Baseline model assuming homogeneous preferences
2. **Latent class logit**: Identifies distinct preference segments
3. **Generalized multinomial logit (G-MNL)**: Accounts for both preference and scale heterogeneity

---

## 8. Pilot Testing Plan

### Phase 1: Qualitative Pre-testing (n = 10-15)
- Think-aloud interviews with target population
- Assess comprehension of attributes, levels, and choice task format
- Identify confusing or unrealistic attribute combinations
- Refine wording and visual layout

### Phase 2: Quantitative Pilot (n = 50-75)
- Administer full survey online
- Assess completion rates, time, dropout patterns
- Evaluate design efficiency with preliminary data
- Update priors for Bayesian D-optimal redesign if needed
- Test data analysis pipeline

### Phase 3: Soft Launch (n = 100)
- First 100 respondents of main survey
- Quality checks: dominance test pass rate, completion time distribution
- Data quality assessment before proceeding to full sample

---

## 9. References

1. Reed Johnson, F., et al. (2013). Constructing experimental designs for discrete-choice experiments: Report of the ISPOR Conjoint Analysis Experimental Design Good Research Practices Task Force. *Value in Health*, 16(1), 3-13.

2. Lancsar, E., & Louviere, J. (2008). Conducting discrete choice experiments to inform healthcare decision making. *PharmacoEconomics*, 26(8), 661-677.

3. de Bekker-Grob, E. W., et al. (2012). Sample size requirements for discrete-choice experiments in healthcare: A practical guide. *Patient*, 8(5), 373-384.

4. Rose, J. M., & Bliemer, M. C. J. (2009). Constructing efficient stated choice experimental designs. *Transport Reviews*, 29(5), 587-617.

5. ChoiceMetrics (2018). *Ngene 1.2 User Manual & Reference Guide*. Sydney, Australia.
