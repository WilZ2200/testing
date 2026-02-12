# Survey Distribution Plan

## Study: Discrete Choice Experiment on Breast Cancer Screening Preferences
### Protocol: DCE-BCS-2026-001 | Version 1.0 | February 2026

---

## 1. Target Population

### 1.1 Primary Population
- **Women aged 40-74 years** residing in the United States
- This age range aligns with current USPSTF breast cancer screening recommendations (biennial mammography for ages 50-74, with individual decision-making for ages 40-49)

### 1.2 Inclusion Criteria
1. Female sex assigned at birth
2. Age 40-74 years at time of survey
3. No current breast cancer diagnosis (active treatment)
4. Able to read and understand English
5. Able to provide informed consent
6. Access to internet (for online administration)

### 1.3 Exclusion Criteria
1. Current breast cancer diagnosis under active treatment
2. Prior bilateral mastectomy
3. Unable to read English at a 6th-grade level
4. Cognitive impairment that prevents informed consent
5. Participation in another breast cancer screening preference study in the past 6 months

### 1.4 Stratification Targets

To ensure a representative sample, we will stratify recruitment by:

| Characteristic | Target Distribution | Source |
|----------------|-------------------|--------|
| Age 40-49 | 30% | US Census proportions |
| Age 50-64 | 45% | US Census proportions |
| Age 65-74 | 25% | US Census proportions |
| White/Caucasian | 58% | US Census |
| Black/African American | 14% | US Census |
| Hispanic/Latina | 18% | US Census |
| Asian/Asian American | 6% | US Census |
| Other races/ethnicities | 4% | US Census |
| Has had prior mammogram | ~65% | NHIS estimates |

---

## 2. Recruitment Channels

### 2.1 Channel 1: Online Research Panel (Primary) -- Prolific

**Target:** 60% of total sample (n = 360)

**Rationale:**
- Higher data quality than MTurk for health-related surveys
- Better demographic diversity and attention
- Participants are naive to choice experiments (less satisficing)
- Built-in screening questionnaire capability
- Fair compensation norms ($8-12/hour)

**Implementation:**
- Create study on Prolific with pre-screening filters (age, sex, country)
- Set target completion time: 20-25 minutes
- Compensation: $5.00 per completion (~$12-15/hour)
- Use Prolific's demographic balancing features
- Estimated recruitment period: 2-3 weeks

### 2.2 Channel 2: Amazon Mechanical Turk (MTurk)

**Target:** 20% of total sample (n = 120)

**Rationale:**
- Large, diverse participant pool
- Fast recruitment
- Lower cost per respondent
- Useful for sensitivity analysis (compare panel effects)

**Implementation:**
- HIT requirements: >98% approval rate, >1,000 approved HITs, US location
- Compensation: $4.00 per completion + $1.00 bonus for quality responses
- Use CloudResearch/TurkPrime for enhanced screening
- Qualification test with basic health literacy questions
- Estimated recruitment period: 1-2 weeks

### 2.3 Channel 3: Clinic-Based Recruitment

**Target:** 15% of total sample (n = 90)

**Rationale:**
- Access to women actively engaged in healthcare
- Higher ecological validity
- Reaches populations underrepresented in online panels
- Validates consistency of online vs. in-person responses

**Implementation:**
- Partner with 3-5 primary care clinics and mammography centers
- Recruit in waiting rooms via tablet-based surveys
- Flyers and provider referrals
- Compensation: $15 gift card per completion
- Research coordinator on-site during recruitment hours
- Estimated recruitment period: 6-8 weeks

### 2.4 Channel 4: Community and Social Media Outreach

**Target:** 5% of total sample (n = 30)

**Rationale:**
- Supplement hard-to-reach demographics
- Reach women who may not be on research panels
- Community organizations provide trusted recruitment channel

**Implementation:**
- Facebook/Instagram targeted ads (women 40-74, health interests)
- Partnerships with community health organizations, breast cancer advocacy groups
- Church/faith-based community outreach
- Compensation: $10 Amazon gift card
- Estimated recruitment period: 4-6 weeks

---

## 3. Sample Size Justification

### 3.1 Rule of Thumb for DCE Sample Size

Using the de Bekker-Grob et al. (2015) formula:

$$n \geq \frac{500 \times c}{t \times a}$$

Where:
- $n$ = minimum number of respondents
- $c$ = largest number of levels for any single attribute = 4 (cost attribute)
- $t$ = number of choice tasks per respondent = 12
- $a$ = number of alternatives per choice set (excluding opt-out) = 2

$$n \geq \frac{500 \times 4}{12 \times 2} = \frac{2000}{24} \approx 84$$

This is the absolute minimum. For a mixed logit model with preference heterogeneity, substantially more observations are needed.

### 3.2 Johnson & Orme (2003) Rule of Thumb

$$n \geq \frac{500 \times c}{t \times a}$$

Same formula, yielding n >= 84 as absolute minimum.

### 3.3 S-estimate from Design

Based on preliminary D-optimal design with Ngene:
- Estimated S-estimate: ~150 respondents for main effects
- For interaction effects and subgroup analysis: ~400-600 respondents

### 3.4 Power Analysis for Subgroup Comparisons

To detect meaningful differences between subgroups (e.g., women with vs. without family history):
- Effect size (Cohen's d): 0.3 (small-medium)
- Power: 80%
- Alpha: 0.05
- Required per group: ~175

### 3.5 Final Target Sample Size

| Component | N |
|-----------|---|
| Minimum for model estimation | 84 |
| Target for mixed logit | 400 |
| Additional for subgroup analysis | +100 |
| Anticipated dropout/data quality exclusions (~15%) | +75 |
| **Final recruitment target** | **600** |

---

## 4. Data Quality Measures

### 4.1 Attention Checks

1. **Dominance test (Choice set 8):** One choice set where Program A is strictly superior on all attributes. Respondents who select Program B or "Neither" are flagged.

2. **Instructed response item:** "To show you are paying attention, please select 'Somewhat difficult' for this question."

3. **Comprehension check:** After the practice choice task, ask: "In the example, which program had the higher sensitivity?" Respondents who answer incorrectly can re-read instructions.

### 4.2 Speeder Detection

- Record time stamps for each choice set
- Flag respondents completing the entire survey in < 8 minutes (estimated minimum reading time)
- Flag respondents completing individual choice sets in < 5 seconds
- Decision rule: Exclude speeders who also fail attention checks

### 4.3 Straightlining Detection

- Identify respondents who choose the same alternative (always A, always B, or always Neither) for all 12 choice sets
- "Always A" or "Always B" in 12/12 sets is suspicious but not automatically exclusionary
- "Always Neither" in 10+ sets suggests disengagement or protest responses

### 4.4 Consistency Test (Test-Retest)

- Choice set 12 is a relabeled repeat of choice set 3
- Compare responses for consistency
- Flag (but do not automatically exclude) inconsistent respondents

### 4.5 Open-ended Response Review

- Review open-ended debriefing responses (D7) for:
  - Evidence of engagement
  - Indicators of confusion
  - "Random" or nonsensical responses
  - Copy-pasted text or bot-generated responses

### 4.6 Data Quality Decision Rules

| Issue | Threshold | Action |
|-------|-----------|--------|
| Failed dominance test | 1 of 1 | Flag; include in sensitivity analysis |
| Speeder (total time) | < 8 min | Exclude |
| Speeder (per choice set) | < 5 sec for 3+ sets | Exclude |
| Complete straightlining | 12/12 same answer | Exclude |
| Failed attention check | 1+ of 2 | Flag; exclude if also speeder |
| Inconsistent test-retest | Different answer | Flag; do not exclude |

**Expected exclusion rate:** 10-15% of raw responses.

---

## 5. Survey Platform and Implementation

### 5.1 Platform: Qualtrics

- Host the survey on Qualtrics with:
  - Randomized block assignment (4 versions)
  - Embedded timing data for each page/question
  - Forced response for choice tasks (prevent skipping)
  - Mobile-responsive design
  - Unique completion codes for panel participants

### 5.2 Technical Requirements

- Browser: Chrome, Firefox, Safari, Edge (latest 2 versions)
- Device: Desktop, tablet, or smartphone (min screen width: 375px)
- Choice task tables must be responsive and readable on mobile
- Progress bar displayed throughout

### 5.3 Survey Flow

1. Eligibility screening (2 questions: age, sex)
2. Informed consent
3. Part A: Screening background (10 questions)
4. Instructions and practice choice task
5. Comprehension check
6. Part B: 12 choice tasks (randomized within block)
7. Part C: Sociodemographics (17 questions)
8. Part D: Debriefing (7 questions)
9. Completion code and thank you page

---

## 6. Timeline

| Phase | Activity | Duration | Start Date |
|-------|----------|----------|------------|
| 1 | Experimental design generation (Ngene) | 2 weeks | March 2026 |
| 2 | Survey programming (Qualtrics) | 2 weeks | March 2026 |
| 3 | IRB submission and approval | 4-6 weeks | April 2026 |
| 4 | Qualitative pre-testing (n=10-15) | 2 weeks | May 2026 |
| 5 | Design revision based on pre-test | 1 week | June 2026 |
| 6 | Quantitative pilot (n=50-75) | 2 weeks | June 2026 |
| 7 | Pilot analysis and design update | 1 week | June 2026 |
| 8 | Full survey launch (soft launch n=100) | 1 week | July 2026 |
| 9 | Full data collection (n=600) | 6-8 weeks | July-August 2026 |
| 10 | Data cleaning and quality checks | 2 weeks | September 2026 |
| 11 | Data analysis | 8 weeks | September-November 2026 |
| 12 | Manuscript preparation | 8 weeks | November 2026-January 2027 |

**Total estimated timeline:** 10-12 months (March 2026 - January 2027)

---

## 7. Budget Estimates

| Item | Unit Cost | Quantity | Total |
|------|-----------|----------|-------|
| **Participant Compensation** | | | |
| Prolific respondents | $5.00 | 400 | $2,000 |
| Prolific platform fee (33%) | | | $660 |
| MTurk respondents | $5.00 | 140 | $700 |
| MTurk platform fee (40%) | | | $280 |
| Clinic-based respondents (gift cards) | $15.00 | 100 | $1,500 |
| Community/social media respondents | $10.00 | 35 | $350 |
| **Software and Platform** | | | |
| Qualtrics (institutional license) | -- | -- | $0* |
| Ngene software license | $595 | 1 | $595 |
| **Recruitment** | | | |
| Facebook/Instagram ad campaign | -- | -- | $500 |
| Printed materials (flyers, consent) | -- | -- | $200 |
| Research coordinator time (clinic) | $25/hr | 80 hrs | $2,000 |
| **Pre-testing** | | | |
| Qualitative interview compensation | $50 | 15 | $750 |
| Pilot study compensation | $5.00 | 75 | $375 |
| | | | |
| **Total estimated budget** | | | **$9,910** |

*Assumes institutional Qualtrics license. If not available, add ~$1,500/year.

---

## 8. Pilot Testing Plan

### 8.1 Phase 1: Cognitive Interviews (n = 10-15)

**Objective:** Assess comprehension, face validity, and usability of the survey instrument.

**Method:**
- Recruit 10-15 women from the target population via convenience sampling
- Conduct 45-60 minute think-aloud interviews (in-person or Zoom)
- Participants complete the survey while verbalizing their thought process
- Interviewer probes on:
  - Understanding of attribute definitions
  - Interpretation of levels (especially sensitivity and false positive rate)
  - Decision-making strategies
  - Layout and readability
  - Time burden and fatigue

**Outcomes:**
- Revised attribute definitions and instructions
- Identification of problematic choice sets
- Improved visual layout

### 8.2 Phase 2: Quantitative Pilot (n = 50-75)

**Objective:** Test survey logistics, timing, and preliminary design performance.

**Method:**
- Administer full survey on Qualtrics via Prolific
- Collect timing data, dropout rates, and response patterns

**Analysis:**
- Median completion time (target: 20-25 minutes)
- Dropout rate by survey section
- Dominance test pass rate (target: >85%)
- Preliminary MNL model estimation to check sign and significance of coefficients
- Design efficiency with real data
- Update Ngene priors if coefficients differ substantially from assumed priors

### 8.3 Phase 3: Soft Launch (n = 100)

**Objective:** Final quality assurance before full data collection.

**Method:**
- First 100 responses from main study
- Implement all data quality checks
- Review open-ended responses

**Go/No-Go Criteria:**
- Completion rate > 80%
- Dominance test pass rate > 85%
- Speeder rate < 15%
- At least 5 of 7 attributes significant in preliminary MNL
- No evidence of systematic technical issues

---

## 9. Ethical Considerations

### 9.1 Informed Consent
- All participants receive a detailed informed consent form before starting
- Consent is documented electronically (checkbox + signature for clinic-based)
- Participants can withdraw at any time without penalty

### 9.2 Privacy and Data Security
- No personally identifiable information (PII) collected in survey
- Prolific/MTurk IDs linked to responses only for compensation; IDs deleted after payment
- Data stored on encrypted, password-protected university server
- Only approved research team members have data access
- Data retained for 7 years per institutional policy

### 9.3 Risk Minimization
- Survey questions do not involve clinical procedures
- Topic may cause mild anxiety about breast cancer; resource links provided at end
- Participants can skip sensitive questions or withdraw
- Debriefing page includes breast cancer screening information resources

### 9.4 Fair Compensation
- All compensation rates meet or exceed platform minimum wage guidelines
- Clinic-based participants receive compensation regardless of survey completion
- No deceptive practices in recruitment or study description
