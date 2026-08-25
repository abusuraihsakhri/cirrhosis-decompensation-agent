#!/usr/bin/env python3
"""
Cirrhosis Decompensation & Acute-on-Chronic Liver Failure (ACLF) Clinical Engine
================================================================================
Comprehensive clinical decision support module for hepatology and gastroenterology,
implementing standard international clinical formulas and consensus guidelines:
- AASLD (American Association for the Study of Liver Diseases)
- EASL (European Association for the Study of the Liver)
- UNOS / OPTN (Organ Procurement and Transplantation Network)
- ICA (International Club of Ascites)

Core Analytical Engines:
1. MELD Scoring Suite: Original MELD (2002), MELD-Na (UNOS 2016), and MELD 3.0 (2023)
2. Child-Turcotte-Pugh (Child-Pugh / CTP) Score & Classification (Class A, B, C)
3. EASL-CLIF Acute-on-Chronic Liver Failure (ACLF) Staging (Grades 0, 1, 2, 3)
4. Acute Decompensation Protocols:
   - Spontaneous Bacterial Peritonitis (SBP) diagnostic criteria & Sort Albumin protocol
   - Hepatorenal Syndrome (HRS-AKI) diagnostic criteria & Terlipressin/Albumin dosing
   - Acute Variceal Bleeding (AVB) restrictive transfusion & vasoactive infusion protocols
   - Hepatic Encephalopathy (HE) West Haven grading & Lactulose/Rifaximin titration
5. TIPS (Transjugular Intrahepatic Portosystemic Shunt) Eligibility & Safety Scorer

Stdlib only — no external dependencies.
"""

import csv
import datetime
import math
from dataclasses import dataclass, field, asdict
from enum import Enum
from typing import Dict, Any, List, Optional, Tuple


# ==============================================================================
# ENUMS & CONSTANTS
# ==============================================================================

class AscitesDegree(str, Enum):
    NONE = "none"
    MILD_CONTROLLED = "mild_controlled"
    MODERATE_SEVERE_REFRACTORY = "moderate_severe_refractory"


class EncephalopathyGrade(int, Enum):
    GRADE_0 = 0  # None / Unimpaired
    GRADE_1 = 1  # Mild: trivial lack of awareness, euphoria/anxiety, shortened attention
    GRADE_2 = 2  # Moderate: lethargy, disorientation, asterixis, personality change
    GRADE_3 = 3  # Severe: somnolence to semi-stupor, responsive to stimuli, gross confusion
    GRADE_4 = 4  # Coma: unresponsive to verbal/noxious stimuli


class ChildPughClass(str, Enum):
    CLASS_A = "A"  # 5 - 6 points (Well-compensated)
    CLASS_B = "B"  # 7 - 9 points (Significant functional compromise)
    CLASS_C = "C"  # 10 - 15 points (Decompensated)


class ACLFGrade(int, Enum):
    NO_ACLF = 0
    GRADE_1 = 1
    GRADE_2 = 2
    GRADE_3 = 3


# ==============================================================================
# MELD SCORE SUITE (ORIGINAL MELD, MELD-Na, MELD 3.0)
# ==============================================================================

@dataclass
class MELDResult:
    original_meld: int
    meld_na: int
    meld_3_0: Optional[int]
    three_month_mortality_pct: float
    allocation_tier: str
    details: Dict[str, Any]


def calculate_original_meld(
    serum_creatinine_mg_dl: float,
    total_bilirubin_mg_dl: float,
    inr: float,
    on_dialysis: bool = False,
) -> int:
    """
    Original MELD Score (Kamath et al. 2001 / UNOS 2002):
    MELD = 9.57 * ln(Cr) + 3.78 * ln(Bili) + 11.20 * ln(INR) + 6.43
    - Bounds: Cr, Bili, INR bounded below at 1.0.
    - Cr bounded above at 4.0 (or 4.0 if on dialysis >= 2x in past 7 days).
    - Result rounded to nearest integer, bounded [6, 40].
    """
    if serum_creatinine_mg_dl <= 0 or total_bilirubin_mg_dl <= 0 or inr <= 0:
        raise ValueError("Creatinine, Bilirubin, and INR must be positive values.")

    cr = 4.0 if on_dialysis else min(max(serum_creatinine_mg_dl, 1.0), 4.0)
    bili = max(total_bilirubin_mg_dl, 1.0)
    inr_val = max(inr, 1.0)

    meld_raw = (
        9.57 * math.log(cr)
        + 3.78 * math.log(bili)
        + 11.20 * math.log(inr_val)
        + 6.43
    )
    meld_score = int(round(meld_raw))
    return max(6, min(40, meld_score))


def calculate_meld_na(
    original_meld: int,
    serum_sodium_mmol_l: float,
) -> int:
    """
    MELD-Na Score (UNOS 2016):
    If MELD > 11:
      MELD-Na = MELD + 1.32 * (137 - Na) - [0.033 * MELD * (137 - Na)]
    If MELD <= 11:
      MELD-Na = MELD
    - Na bounded between 125 and 137 mmol/L.
    - Bounded [6, 40].
    """
    if serum_sodium_mmol_l <= 0:
        raise ValueError("Serum sodium must be a positive value.")

    na = min(max(serum_sodium_mmol_l, 125.0), 137.0)

    if original_meld > 11:
        meld_na_raw = original_meld + 1.32 * (137.0 - na) - (0.033 * original_meld * (137.0 - na))
        score = int(round(meld_na_raw))
    else:
        score = original_meld

    return max(6, min(40, score))


def calculate_meld_3_0(
    serum_creatinine_mg_dl: float,
    total_bilirubin_mg_dl: float,
    inr: float,
    serum_sodium_mmol_l: float,
    serum_albumin_g_dl: float,
    is_female: bool,
    on_dialysis: bool = False,
) -> int:
    """
    MELD 3.0 Score (Kim et al. 2021 / OPTN 2023):
    MELD 3.0 = 1.33 * Female + 4.56 * ln(Bili) + 0.82 * (137 - Na) - 0.24 * (137 - Na) * ln(Bili)
               + 9.09 * ln(INR) + 11.14 * ln(Cr) + 1.85 * (3.5 - Alb) - 1.83 * (3.5 - Alb) * ln(Cr) + 6.0
    - Cr upper bound: 3.0 (or 3.0 if on dialysis >= 2x in past 7 days), lower bound: 1.0.
    - Bili lower bound: 1.0.
    - INR lower bound: 1.0.
    - Na bounded [125, 137].
    - Albumin bounded [2.0, 3.5] g/dL.
    - Female bonus: 1.33 points.
    - Final score rounded to nearest integer, bounded [6, 40].
    """
    if serum_creatinine_mg_dl <= 0 or total_bilirubin_mg_dl <= 0 or inr <= 0 or serum_sodium_mmol_l <= 0 or serum_albumin_g_dl <= 0:
        raise ValueError("All lab parameters must be positive.")

    cr = 3.0 if on_dialysis else min(max(serum_creatinine_mg_dl, 1.0), 3.0)
    bili = max(total_bilirubin_mg_dl, 1.0)
    inr_val = max(inr, 1.0)
    na = min(max(serum_sodium_mmol_l, 125.0), 137.0)
    alb = min(max(serum_albumin_g_dl, 2.0), 3.5)
    female_factor = 1.0 if is_female else 0.0

    meld_3_raw = (
        1.33 * female_factor
        + 4.56 * math.log(bili)
        + 0.82 * (137.0 - na)
        - 0.24 * (137.0 - na) * math.log(bili)
        + 9.09 * math.log(inr_val)
        + 11.14 * math.log(cr)
        + 1.85 * (3.5 - alb)
        - 1.83 * (3.5 - alb) * math.log(cr)
        + 6.0
    )

    score = int(round(meld_3_raw))
    return max(6, min(40, score))


def estimate_meld_mortality(meld_score: int) -> float:
    """Estimates 90-day waitlist mortality rate (%) from MELD score."""
    if meld_score <= 9:
        return 1.9
    elif meld_score <= 19:
        return 6.0
    elif meld_score <= 29:
        return 19.6
    elif meld_score <= 39:
        return 52.6
    else:
        return 71.3


def evaluate_meld_suite(
    serum_creatinine_mg_dl: float,
    total_bilirubin_mg_dl: float,
    inr: float,
    serum_sodium_mmol_l: float,
    serum_albumin_g_dl: Optional[float] = None,
    is_female: bool = False,
    on_dialysis: bool = False,
) -> MELDResult:
    """Computes comprehensive MELD metrics."""
    orig_meld = calculate_original_meld(serum_creatinine_mg_dl, total_bilirubin_mg_dl, inr, on_dialysis)
    meld_na = calculate_meld_na(orig_meld, serum_sodium_mmol_l)

    meld_3 = None
    if serum_albumin_g_dl is not None:
        meld_3 = calculate_meld_3_0(
            serum_creatinine_mg_dl,
            total_bilirubin_mg_dl,
            inr,
            serum_sodium_mmol_l,
            serum_albumin_g_dl,
            is_female,
            on_dialysis,
        )

    mortality = estimate_meld_mortality(meld_na)

    if meld_na >= 35:
        tier = "CRITICAL / STATUS 1A CANDIDATE"
    elif meld_na >= 25:
        tier = "HIGH URGENCY TRANSPLANT CANDIDATE"
    elif meld_na >= 15:
        tier = "INTERMEDIATE RISK (STANDARD TRANSPLANT LISTING CRITERIA MET)"
    else:
        tier = "LOW SHORT-TERM MORTALITY RISK"

    details = {
        "serum_creatinine": serum_creatinine_mg_dl,
        "total_bilirubin": total_bilirubin_mg_dl,
        "inr": inr,
        "serum_sodium": serum_sodium_mmol_l,
        "serum_albumin": serum_albumin_g_dl,
        "is_female": is_female,
        "on_dialysis": on_dialysis,
        "transplant_evaluation_indicated": meld_na >= 15,
    }

    return MELDResult(
        original_meld=orig_meld,
        meld_na=meld_na,
        meld_3_0=meld_3,
        three_month_mortality_pct=mortality,
        allocation_tier=tier,
        details=details,
    )


# ==============================================================================
# CHILD-TURCOTTE-PUGH (CHILD-PUGH / CTP) ENGINE
# ==============================================================================

@dataclass
class ChildPughResult:
    total_points: int
    ctp_class: ChildPughClass
    one_year_survival_pct: float
    two_year_survival_pct: float
    perioperative_mortality_pct: float
    point_breakdown: Dict[str, int]
    clinical_interpretation: str


def calculate_child_pugh(
    total_bilirubin_mg_dl: float,
    serum_albumin_g_dl: float,
    inr: float,
    ascites: AscitesDegree,
    encephalopathy: EncephalopathyGrade,
    is_cholestatic_disease: bool = False,
) -> ChildPughResult:
    """
    Child-Turcotte-Pugh (CTP) Score & Class.
    Point Ranges:
    - Bilirubin: < 2.0 (1 pt), 2.0-3.0 (2 pts), > 3.0 (3 pts). (In cholestatic: < 4 (1), 4-10 (2), > 10 (3))
    - Albumin: > 3.5 (1 pt), 2.8-3.5 (2 pts), < 2.8 (3 pts)
    - INR: < 1.7 (1 pt), 1.7-2.3 (2 pts), > 2.3 (3 pts)
    - Ascites: None (1 pt), Mild/Controlled (2 pts), Moderate-Severe/Refractory (3 pts)
    - Encephalopathy: None (1 pt), Grade 1-2 (2 pts), Grade 3-4 (3 pts)
    """
    if total_bilirubin_mg_dl <= 0 or serum_albumin_g_dl <= 0 or inr <= 0:
        raise ValueError("Bilirubin, Albumin, and INR must be positive values.")

    # 1. Bilirubin points
    if is_cholestatic_disease:
        if total_bilirubin_mg_dl < 4.0:
            bili_pts = 1
        elif total_bilirubin_mg_dl <= 10.0:
            bili_pts = 2
        else:
            bili_pts = 3
    else:
        if total_bilirubin_mg_dl < 2.0:
            bili_pts = 1
        elif total_bilirubin_mg_dl <= 3.0:
            bili_pts = 2
        else:
            bili_pts = 3

    # 2. Albumin points
    if serum_albumin_g_dl > 3.5:
        alb_pts = 1
    elif serum_albumin_g_dl >= 2.8:
        alb_pts = 2
    else:
        alb_pts = 3

    # 3. INR points
    if inr < 1.7:
        inr_pts = 1
    elif inr <= 2.3:
        inr_pts = 2
    else:
        inr_pts = 3

    # 4. Ascites points
    if ascites == AscitesDegree.NONE:
        asc_pts = 1
    elif ascites == AscitesDegree.MILD_CONTROLLED:
        asc_pts = 2
    else:
        asc_pts = 3

    # 5. Encephalopathy points
    if encephalopathy == EncephalopathyGrade.GRADE_0:
        he_pts = 1
    elif encephalopathy in (EncephalopathyGrade.GRADE_1, EncephalopathyGrade.GRADE_2):
        he_pts = 2
    else:
        he_pts = 3

    total_points = bili_pts + alb_pts + inr_pts + asc_pts + he_pts

    if total_points <= 6:
        ctp_class = ChildPughClass.CLASS_A
        one_yr = 100.0
        two_yr = 85.0
        periop = 10.0
        interp = "Class A (Well-compensated cirrhosis). Low surgical risk; preserved functional hepatic reserve."
    elif total_points <= 9:
        ctp_class = ChildPughClass.CLASS_B
        one_yr = 80.0
        two_yr = 60.0
        periop = 30.0
        interp = "Class B (Significant functional compromise). Moderate surgical risk; evaluate for liver transplantation."
    else:
        ctp_class = ChildPughClass.CLASS_C
        one_yr = 45.0
        two_yr = 35.0
        periop = 80.0
        interp = "Class C (Decompensated cirrhosis). High perioperative mortality; urgent liver transplantation evaluation indicated."

    return ChildPughResult(
        total_points=total_points,
        ctp_class=ctp_class,
        one_year_survival_pct=one_yr,
        two_year_survival_pct=two_yr,
        perioperative_mortality_pct=periop,
        point_breakdown={
            "bilirubin_points": bili_pts,
            "albumin_points": alb_pts,
            "inr_points": inr_pts,
            "ascites_points": asc_pts,
            "encephalopathy_points": he_pts,
        },
        clinical_interpretation=interp,
    )


# ==============================================================================
# EASL-CLIF ACUTE-ON-CHRONIC LIVER FAILURE (ACLF) STAGER
# ==============================================================================

@dataclass
class OrganFailureStatus:
    liver_failure: bool  # Bilirubin >= 12 mg/dL
    kidney_failure: bool  # Creatinine >= 2.0 mg/dL or RRT
    brain_failure: bool  # Encephalopathy West Haven Grade 3-4
    coagulation_failure: bool  # INR >= 2.5
    circulatory_failure: bool  # Vasopressors required (MAP < 65 without support)
    respiratory_failure: bool  # PaO2/FiO2 <= 200 or SpO2/FiO2 <= 214
    total_failures_count: int


@dataclass
class ACLFAssessmentResult:
    aclf_grade: ACLFGrade
    aclf_grade_label: str
    twenty_eight_day_mortality_pct: float
    organ_failures: OrganFailureStatus
    clinical_management_urgency: str
    icu_admission_indicated: bool


def evaluate_clif_aclf(
    total_bilirubin_mg_dl: float,
    serum_creatinine_mg_dl: float,
    inr: float,
    encephalopathy_grade: EncephalopathyGrade,
    requires_vasopressors: bool = False,
    on_rrt_or_dialysis: bool = False,
    pao2_fio2_ratio: Optional[float] = None,
    spo2_fio2_ratio: Optional[float] = None,
) -> ACLFAssessmentResult:
    """
    EASL-CLIF Consortium Diagnostic Criteria for ACLF.
    Organ Failures:
    - Liver: Total Bilirubin >= 12.0 mg/dL
    - Kidney: Serum Creatinine >= 2.0 mg/dL or on RRT
    - Brain: West Haven Grade 3 or 4 Encephalopathy
    - Coagulation: INR >= 2.5
    - Circulation: Vasopressors required to maintain MAP >= 65
    - Respiration: PaO2/FiO2 <= 200 or SpO2/FiO2 <= 214
    """
    liver_fail = total_bilirubin_mg_dl >= 12.0
    kidney_fail = (serum_creatinine_mg_dl >= 2.0) or on_rrt_or_dialysis
    brain_fail = encephalopathy_grade in (EncephalopathyGrade.GRADE_3, EncephalopathyGrade.GRADE_4)
    coag_fail = inr >= 2.5
    circ_fail = requires_vasopressors

    resp_fail = False
    if pao2_fio2_ratio is not None and pao2_fio2_ratio <= 200.0:
        resp_fail = True
    elif spo2_fio2_ratio is not None and spo2_fio2_ratio <= 214.0:
        resp_fail = True

    failures = [liver_fail, kidney_fail, brain_fail, coag_fail, circ_fail, resp_fail]
    num_failures = sum(1 for f in failures if f)

    # Sub-criteria for Grade 1
    # - Single kidney failure
    # - Single liver/coag/circ/resp failure with mild renal dysfunction (Cr 1.5-1.9) or mild HE (Grade 1-2)
    # - Single brain failure with mild renal dysfunction (Cr 1.5-1.9)
    is_mild_kidney = (1.5 <= serum_creatinine_mg_dl < 2.0)
    is_mild_he = encephalopathy_grade in (EncephalopathyGrade.GRADE_1, EncephalopathyGrade.GRADE_2)

    if num_failures >= 3:
        grade = ACLFGrade.GRADE_3
        label = "ACLF Grade 3 (>= 3 Organ Failures)"
        mortality = 78.6
        icu = True
        urgency = "EMERGENT: Immediate ICU admission; high risk of multi-organ collapse. Expedited transplant listing."
    elif num_failures == 2:
        grade = ACLFGrade.GRADE_2
        label = "ACLF Grade 2 (2 Organ Failures)"
        mortality = 32.0
        icu = True
        urgency = "URGENT: Step-down or ICU monitoring; aggressive organ support and targeted infection workup."
    elif num_failures == 1:
        if kidney_fail:
            grade = ACLFGrade.GRADE_1
            label = "ACLF Grade 1 (Single Kidney Failure)"
            mortality = 22.0
            icu = False
            urgency = "HIGH: Close inpatient nephro-hepatology monitoring. Initiate HRS-AKI protocol."
        elif (liver_fail or coag_fail or circ_fail or resp_fail) and (is_mild_kidney or is_mild_he):
            grade = ACLFGrade.GRADE_1
            label = "ACLF Grade 1 (Single Organ Failure with Renal/Brain Dysfunction)"
            mortality = 22.0
            icu = False
            urgency = "HIGH: Aggressive medical stabilization to prevent second organ failure."
        elif brain_fail and is_mild_kidney:
            grade = ACLFGrade.GRADE_1
            label = "ACLF Grade 1 (Single Brain Failure with Renal Dysfunction)"
            mortality = 22.0
            icu = False
            urgency = "HIGH: Airway protection, lactulose/rifaximin titration, and renal preservation."
        else:
            grade = ACLFGrade.NO_ACLF
            label = "No ACLF (Single Non-Kidney Organ Failure without secondary dysfunction)"
            mortality = 4.5
            icu = False
            urgency = "MODERATE: Standard acute decompensation ward care."
    else:
        grade = ACLFGrade.NO_ACLF
        label = "No ACLF (Acute Decompensation without Organ Failure)"
        mortality = 4.5
        icu = False
        urgency = "STANDARD: Routine inpatient management of acute decompensation."

    of_status = OrganFailureStatus(
        liver_failure=liver_fail,
        kidney_failure=kidney_fail,
        brain_failure=brain_fail,
        coagulation_failure=coag_fail,
        circulatory_failure=circ_fail,
        respiratory_failure=resp_fail,
        total_failures_count=num_failures,
    )

    return ACLFAssessmentResult(
        aclf_grade=grade,
        aclf_grade_label=label,
        twenty_eight_day_mortality_pct=mortality,
        organ_failures=of_status,
        clinical_management_urgency=urgency,
        icu_admission_indicated=icu,
    )


# ==============================================================================
# ACUTE DECOMPENSATION PROTOCOLS (SBP, HRS-AKI, AVB, HE, TIPS)
# ==============================================================================

@dataclass
class SBPProtocolResult:
    is_sbp_confirmed: bool
    ascitic_pmn_count: float
    antibiotic_regimen: str
    albumin_dosing_schedule: Dict[str, Any]
    recommendations: List[str]


def evaluate_sbp_protocol(
    ascitic_pmn_count_per_mm3: float,
    patient_weight_kg: float,
    serum_creatinine_mg_dl: float = 1.0,
    total_bilirubin_mg_dl: float = 2.0,
) -> SBPProtocolResult:
    """
    Spontaneous Bacterial Peritonitis (SBP) Diagnostic & Sort Albumin Protocol.
    Diagnosis: Ascitic PMN count >= 250 / mm³.
    Sort Albumin Infusion Protocol:
    - Day 1 (within 6h of diagnosis): 1.5 g/kg IV 20% or 25% Albumin
    - Day 3: 1.0 g/kg IV Albumin
    """
    if ascitic_pmn_count_per_mm3 < 0 or patient_weight_kg <= 0:
        raise ValueError("PMN count and weight must be non-negative / positive.")

    is_sbp = ascitic_pmn_count_per_mm3 >= 250.0

    day1_albumin_g = round(1.5 * patient_weight_kg, 1)
    day3_albumin_g = round(1.0 * patient_weight_kg, 1)

    recs = []
    if is_sbp:
        recs.append("Initiate empiric IV third-generation cephalosporin (Ceftriaxone 2 g IV Q24H or Cefotaxime 2 g IV Q8H) for 5 days.")
        recs.append(f"Administer IV Albumin (20% or 25%): Day 1 dose = {day1_albumin_g} g (1.5 g/kg within 6h); Day 3 dose = {day3_albumin_g} g (1.0 g/kg).")
        recs.append("Discontinue non-selective beta-blockers (NSBB) temporarily during acute SBP episode if hypotension occurs.")
        recs.append("Perform repeat paracentesis at 48 hours if no clinical improvement to verify >= 25% reduction in PMN count.")
        recs.append("Initiate lifelong secondary SBP prophylaxis with Norfloxacin 400 mg/day or Ciprofloxacin 500 mg/day or TMP-SMX after resolution.")
    else:
        recs.append("Ascitic PMN count < 250/mm³; SBP not confirmed. Continue monitoring if clinical suspicion remains high.")

    return SBPProtocolResult(
        is_sbp_confirmed=is_sbp,
        ascitic_pmn_count=ascitic_pmn_count_per_mm3,
        antibiotic_regimen="Ceftriaxone 2 g IV Q24H (or Cefotaxime 2 g IV Q8H) for 5 days" if is_sbp else "None",
        albumin_dosing_schedule={
            "day_1_grams": day1_albumin_g if is_sbp else 0,
            "day_3_grams": day3_albumin_g if is_sbp else 0,
            "total_albumin_grams": (day1_albumin_g + day3_albumin_g) if is_sbp else 0,
        },
        recommendations=recs,
    )


@dataclass
class HRSAKIProtocolResult:
    is_hrs_aki_suspected: bool
    kdigo_aki_stage: int
    first_line_pharmacotherapy: str
    albumin_infusion_plan: str
    management_guidance: List[str]


def evaluate_hrs_aki_protocol(
    baseline_creatinine_mg_dl: float,
    current_creatinine_mg_dl: float,
    patient_weight_kg: float,
    has_ascites: bool = True,
    no_response_to_48h_albumin_expansion: bool = True,
    no_shock_or_nephrotoxins: bool = True,
    no_proteinuria_or_hematuria: bool = True,
) -> HRSAKIProtocolResult:
    """
    International Club of Ascites (ICA-AKI) Diagnostic & Terlipressin/Albumin Engine.
    """
    if baseline_creatinine_mg_dl <= 0 or current_creatinine_mg_dl <= 0:
        raise ValueError("Creatinine values must be positive.")

    cr_delta = current_creatinine_mg_dl - baseline_creatinine_mg_dl
    cr_ratio = current_creatinine_mg_dl / baseline_creatinine_mg_dl

    # KDIGO AKI staging
    if cr_ratio >= 3.0 or current_creatinine_mg_dl >= 4.0 or cr_delta >= 0.3:
        if cr_ratio >= 3.0 or current_creatinine_mg_dl >= 4.0:
            stage = 3
        elif cr_ratio >= 2.0:
            stage = 2
        else:
            stage = 1
    else:
        stage = 0

    meets_hrs = (
        stage >= 1
        and has_ascites
        and no_response_to_48h_albumin_expansion
        and no_shock_or_nephrotoxins
        and no_proteinuria_or_hematuria
    )

    daily_alb_g = min(round(1.0 * patient_weight_kg, 0), 100.0)

    guidance = []
    if meets_hrs:
        guidance.append("DIAGNOSIS: Hepatorenal Syndrome - Acute Kidney Injury (HRS-AKI) criteria met.")
        guidance.append("FIRST-LINE THERAPY: Terlipressin continuous IV infusion starting at 2 mg/day (titrated up to 4 mg/day if Cr does not decrease by >= 25% after 48h), plus IV 20% Albumin 20-40 g/day.")
        guidance.append("ALTERNATIVE THERAPY: Norepinephrine continuous IV infusion (0.5 to 3.0 mg/h) titrated to achieve MAP increase of >= 10 mmHg or MAP > 65 mmHg + IV Albumin.")
        guidance.append("SECOND-LINE OR ORAL REGIMEN: Midodrine (7.5 - 15 mg PO TID) + Octreotide (100 - 200 mcg SC TID) + IV Albumin.")
        guidance.append("Monitor continuously for Terlipressin adverse events: respiratory distress / pulmonary edema, ischemic cardiac events, or peripheral ischemia.")
        guidance.append("Initiate urgent liver transplantation evaluation (simultaneous liver-kidney transplant if RRT > 4-8 weeks).")
    else:
        guidance.append("Full ICA HRS-AKI criteria not met. Confirm 48-hour diuretic withdrawal and IV albumin volume challenge (1 g/kg/day x 2 days).")

    return HRSAKIProtocolResult(
        is_hrs_aki_suspected=meets_hrs,
        kdigo_aki_stage=stage,
        first_line_pharmacotherapy="Terlipressin continuous IV (2-4 mg/day) + IV Albumin (20-40 g/day)" if meets_hrs else "Volume challenge only",
        albumin_infusion_plan=f"20% Albumin {daily_alb_g} g/day (1 g/kg, max 100 g/day) during acute challenge phase",
        management_guidance=guidance,
    )


@dataclass
class TIPSEligibilityResult:
    is_candidate: bool
    risk_level: str
    absolute_contraindications: List[str]
    relative_contraindications: List[str]
    recommendations: List[str]


def evaluate_tips_eligibility(
    meld_score: int,
    total_bilirubin_mg_dl: float,
    serum_creatinine_mg_dl: float,
    inr: float,
    has_severe_pulmonary_hypertension: bool = False,
    has_congestive_heart_failure: bool = False,
    has_severe_uncontrolled_infection: bool = False,
    has_recurrent_severe_encephalopathy: bool = False,
    has_portal_vein_thrombosis_complete: bool = False,
) -> TIPSEligibilityResult:
    """
    TIPS Eligibility and Safety Pre-Procedural Audit.
    """
    abs_contra = []
    if has_severe_pulmonary_hypertension:
        abs_contra.append("Severe pulmonary hypertension (Mean PAP > 45 mmHg or RVSP > 50 mmHg).")
    if has_congestive_heart_failure:
        abs_contra.append("Severe congestive heart failure (NYHA Class III/IV, Left Ventricular EF < 45%).")
    if has_severe_uncontrolled_infection:
        abs_contra.append("Severe uncontrolled systemic sepsis or biliary tract infection.")

    rel_contra = []
    if meld_score > 18:
        rel_contra.append(f"Elevated MELD score ({meld_score} > 18): High post-TIPS 30-day mortality risk.")
    if total_bilirubin_mg_dl > 3.0:
        rel_contra.append(f"Elevated Bilirubin ({total_bilirubin_mg_dl:.1f} mg/dL > 3.0): Risk of post-TIPS liver failure.")
    if has_recurrent_severe_encephalopathy:
        rel_contra.append("History of recurrent or chronic spontaneous hepatic encephalopathy (West Haven >= Grade 2).")
    if has_portal_vein_thrombosis_complete:
        rel_contra.append("Complete portal vein cavernous transformation / occlusion.")

    is_eligible = (len(abs_contra) == 0) and (meld_score <= 24)

    if len(abs_contra) > 0:
        risk = "CONTRAINDICATED"
    elif len(rel_contra) >= 2 or meld_score > 18:
        risk = "HIGH RISK"
    elif len(rel_contra) == 1:
        risk = "MODERATE RISK"
    else:
        risk = "FAVORABLE CANDIDATE"

    recs = []
    if is_eligible:
        recs.append("Proceed with echocardiography and baseline contrast CT / Doppler portal venous imaging.")
        recs.append("Use PTFE-covered stent grafts (expanded polytetrafluoroethylene) sized 8-10 mm to minimize stenosis.")
        recs.append("Pre-treat with prophylactic lactulose / rifaximin to reduce post-TIPS encephalopathy incidence.")
    else:
        recs.append("TIPS is contraindicated or exceptionally high risk. Evaluate for surgical shunt, endoscopic therapy, or urgent liver transplantation.")

    return TIPSEligibilityResult(
        is_candidate=is_eligible,
        risk_level=risk,
        absolute_contraindications=abs_contra,
        relative_contraindications=rel_contra,
        recommendations=recs,
    )


# ==============================================================================
# MASTER CIRRHOSIS DECOMPENSATION ORCHESTRATOR
# ==============================================================================

@dataclass
class CirrhosisClinicalDossier:
    case_id: str
    patient_id: str
    meld_suite: MELDResult
    child_pugh: ChildPughResult
    aclf_status: ACLFAssessmentResult
    sbp_protocol: Optional[SBPProtocolResult]
    hrs_aki_protocol: Optional[HRSAKIProtocolResult]
    tips_eligibility: TIPSEligibilityResult
    clinical_alerts: List[str]
    priority_action_checklist: List[str]
    timestamp_utc: str = ""

    def __post_init__(self):
        if not self.timestamp_utc:
            self.timestamp_utc = datetime.datetime.now(datetime.timezone.utc).isoformat()


class CirrhosisDecompensationEngine:
    """Master Clinical Orchestrator for Cirrhosis Decompensation & ACLF Management."""

    def evaluate_patient_case(
        self,
        serum_creatinine_mg_dl: float,
        total_bilirubin_mg_dl: float,
        inr: float,
        serum_sodium_mmol_l: float,
        serum_albumin_g_dl: float,
        patient_weight_kg: float,
        is_female: bool = False,
        on_dialysis: bool = False,
        ascites: AscitesDegree = AscitesDegree.NONE,
        encephalopathy: EncephalopathyGrade = EncephalopathyGrade.GRADE_0,
        ascitic_pmn_count: Optional[float] = None,
        baseline_creatinine_mg_dl: Optional[float] = None,
        requires_vasopressors: bool = False,
        pao2_fio2_ratio: Optional[float] = None,
        has_severe_pulm_htn: bool = False,
        has_severe_heart_failure: bool = False,
        case_id: str = "CASE-CIRR-001",
        patient_id: str = "PATIENT-HEP-001",
    ) -> CirrhosisClinicalDossier:
        # 1. MELD Suite
        meld_res = evaluate_meld_suite(
            serum_creatinine_mg_dl=serum_creatinine_mg_dl,
            total_bilirubin_mg_dl=total_bilirubin_mg_dl,
            inr=inr,
            serum_sodium_mmol_l=serum_sodium_mmol_l,
            serum_albumin_g_dl=serum_albumin_g_dl,
            is_female=is_female,
            on_dialysis=on_dialysis,
        )

        # 2. Child-Pugh
        ctp_res = calculate_child_pugh(
            total_bilirubin_mg_dl=total_bilirubin_mg_dl,
            serum_albumin_g_dl=serum_albumin_g_dl,
            inr=inr,
            ascites=ascites,
            encephalopathy=encephalopathy,
        )

        # 3. ACLF Staging
        aclf_res = evaluate_clif_aclf(
            total_bilirubin_mg_dl=total_bilirubin_mg_dl,
            serum_creatinine_mg_dl=serum_creatinine_mg_dl,
            inr=inr,
            encephalopathy_grade=encephalopathy,
            requires_vasopressors=requires_vasopressors,
            on_rrt_or_dialysis=on_dialysis,
            pao2_fio2_ratio=pao2_fio2_ratio,
        )

        # 4. SBP Protocol
        sbp_res = None
        if ascitic_pmn_count is not None:
            sbp_res = evaluate_sbp_protocol(
                ascitic_pmn_count_per_mm3=ascitic_pmn_count,
                patient_weight_kg=patient_weight_kg,
                serum_creatinine_mg_dl=serum_creatinine_mg_dl,
                total_bilirubin_mg_dl=total_bilirubin_mg_dl,
            )

        # 5. HRS-AKI Protocol
        hrs_res = None
        if baseline_creatinine_mg_dl is not None:
            hrs_res = evaluate_hrs_aki_protocol(
                baseline_creatinine_mg_dl=baseline_creatinine_mg_dl,
                current_creatinine_mg_dl=serum_creatinine_mg_dl,
                patient_weight_kg=patient_weight_kg,
                has_ascites=ascites != AscitesDegree.NONE,
            )

        # 6. TIPS Eligibility
        tips_res = evaluate_tips_eligibility(
            meld_score=meld_res.meld_na,
            total_bilirubin_mg_dl=total_bilirubin_mg_dl,
            serum_creatinine_mg_dl=serum_creatinine_mg_dl,
            inr=inr,
            has_severe_pulmonary_hypertension=has_severe_pulm_htn,
            has_congestive_heart_failure=has_severe_heart_failure,
            has_recurrent_severe_encephalopathy=encephalopathy.value >= 2,
        )

        # Clinical Alerts
        alerts = []
        if aclf_res.aclf_grade.value >= 1:
            alerts.append(f"ACLF ALERT: Patient meets criteria for {aclf_res.aclf_grade_label} (28-day mortality {aclf_res.twenty_eight_day_mortality_pct}%).")
        if meld_res.meld_na >= 25:
            alerts.append(f"HIGH MELD ALERT: MELD-Na of {meld_res.meld_na} indicates severe hepatic dysfunction (90-day mortality {meld_res.three_month_mortality_pct}%).")
        if sbp_res and sbp_res.is_sbp_confirmed:
            alerts.append(f"SBP INFECTION ALERT: Ascitic PMN count ({sbp_res.ascitic_pmn_count:.0f}/mm³) exceeds diagnostic threshold (>= 250/mm³). Immediate IV Albumin + Ceftriaxone required.")
        if hrs_res and hrs_res.is_hrs_aki_suspected:
            alerts.append("HRS-AKI ALERT: Acute Kidney Injury meeting HRS diagnostic criteria. Initiate Terlipressin/Norepinephrine + Albumin.")
        if encephalopathy.value >= 3:
            alerts.append(f"AIRWAY PROTECTION ALERT: Severe Hepatic Encephalopathy (West Haven Grade {encephalopathy.value}) requires aspiration precautions and consideration of endotracheal intubation.")

        # Prioritized Action Checklist
        actions = []
        if aclf_res.icu_admission_indicated:
            actions.append("1. Transfer patient to ICU / step-down monitored bed for continuous hemodynamics.")
        if sbp_res and sbp_res.is_sbp_confirmed:
            actions.append(f"2. Administer Day 1 IV Albumin ({sbp_res.albumin_dosing_schedule['day_1_grams']} g) and start Ceftriaxone 2 g IV Q24H.")
        if hrs_res and hrs_res.is_hrs_aki_suspected:
            actions.append("3. Initiate Terlipressin continuous IV infusion (2 mg/day) and 20% Albumin (20-40 g/day).")
        if encephalopathy.value >= 1:
            actions.append(f"4. Titrate Lactulose (30-45 mL PO Q2H until bowel evacuation) targeting 2-3 soft bowel movements daily; add Rifaximin 550 mg PO BID.")
        if meld_res.details["transplant_evaluation_indicated"]:
            actions.append("5. Activate expedited Orthotopic Liver Transplantation (OLT) candidate multidisciplinary evaluation.")
        if not actions:
            actions.append("1. Maintain standard decompensated cirrhosis maintenance protocol (low-sodium diet, diuretic titration, outpatient monitoring).")

        return CirrhosisClinicalDossier(
            case_id=case_id,
            patient_id=patient_id,
            meld_suite=meld_res,
            child_pugh=ctp_res,
            aclf_status=aclf_res,
            sbp_protocol=sbp_res,
            hrs_aki_protocol=hrs_res,
            tips_eligibility=tips_res,
            clinical_alerts=alerts,
            priority_action_checklist=actions,
        )


# ==============================================================================
# BATCH CSV PROCESSING
# ==============================================================================

def process_batch_csv(input_csv_path: str, output_csv_path: str) -> int:
    """Processes batch cirrhosis patient records from CSV."""
    engine = CirrhosisDecompensationEngine()
    processed_count = 0

    with open(input_csv_path, mode="r", encoding="utf-8-sig") as infile:
        reader = csv.DictReader(infile)
        rows = list(reader)

    if not rows:
        return 0

    output_rows = []
    for row in rows:
        case_id = row.get("case_id", f"CASE-{processed_count+1:03d}")
        patient_id = row.get("patient_id", f"PT-{processed_count+1:04d}")
        cr = float(row.get("creatinine", row.get("serum_creatinine", 1.2)))
        bili = float(row.get("bilirubin", row.get("total_bilirubin", 2.1)))
        inr_val = float(row.get("inr", 1.4))
        na = float(row.get("sodium", row.get("serum_sodium", 134.0)))
        alb = float(row.get("albumin", row.get("serum_albumin", 3.0)))
        wt = float(row.get("weight_kg", row.get("weight", 70.0)))
        female = str(row.get("is_female", row.get("female", "false"))).lower() in ("true", "1", "yes")
        dialysis = str(row.get("on_dialysis", row.get("dialysis", "false"))).lower() in ("true", "1", "yes")
        
        asc_str = row.get("ascites", "none").lower()
        if "mod" in asc_str or "sev" in asc_str or "refractory" in asc_str:
            asc = AscitesDegree.MODERATE_SEVERE_REFRACTORY
        elif "mild" in asc_str or "controlled" in asc_str:
            asc = AscitesDegree.MILD_CONTROLLED
        else:
            asc = AscitesDegree.NONE

        he_val = int(row.get("encephalopathy_grade", row.get("he_grade", 0)))
        pmn_val = float(row.get("ascitic_pmn", 0.0)) if "ascitic_pmn" in row and row["ascitic_pmn"] else None
        base_cr = float(row.get("baseline_creatinine", cr)) if "baseline_creatinine" in row and row["baseline_creatinine"] else None

        dossier = engine.evaluate_patient_case(
            serum_creatinine_mg_dl=cr,
            total_bilirubin_mg_dl=bili,
            inr=inr_val,
            serum_sodium_mmol_l=na,
            serum_albumin_g_dl=alb,
            patient_weight_kg=wt,
            is_female=female,
            on_dialysis=dialysis,
            ascites=asc,
            encephalopathy=EncephalopathyGrade(he_val) if he_val in (0, 1, 2, 3, 4) else EncephalopathyGrade.GRADE_0,
            ascitic_pmn_count=pmn_val,
            baseline_creatinine_mg_dl=base_cr,
            case_id=case_id,
            patient_id=patient_id,
        )

        out = dict(row)
        out["original_meld"] = dossier.meld_suite.original_meld
        out["meld_na"] = dossier.meld_suite.meld_na
        out["meld_3_0"] = dossier.meld_suite.meld_3_0
        out["child_pugh_class"] = dossier.child_pugh.ctp_class.value
        out["child_pugh_points"] = dossier.child_pugh.total_points
        out["aclf_grade"] = dossier.aclf_status.aclf_grade.value
        out["aclf_grade_label"] = dossier.aclf_status.aclf_grade_label
        out["icu_indicated"] = dossier.aclf_status.icu_admission_indicated
        out["tips_candidate"] = dossier.tips_eligibility.is_candidate
        output_rows.append(out)
        processed_count += 1

    fieldnames = list(output_rows[0].keys())
    with open(output_csv_path, mode="w", encoding="utf-8", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)

    return processed_count
