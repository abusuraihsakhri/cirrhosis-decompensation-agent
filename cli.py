#!/usr/bin/env python3
"""
Command-Line Interface for Cirrhosis Decompensation & ACLF Clinical Decision Support Engine.

Usage:
    python cli.py evaluate --cr 1.8 --bili 3.5 --inr 1.9 --na 132 --alb 2.6 --weight 75 --female --ascites moderate --he 2 --pmn 320
    python cli.py meld --cr 2.1 --bili 4.2 --inr 2.0 --na 130 --alb 2.7 --female
    python cli.py child-pugh --bili 3.2 --alb 2.5 --inr 2.4 --ascites moderate --he 2
    python cli.py aclf --cr 2.4 --bili 14.0 --inr 2.6 --he 3 --vasopressors
    python cli.py sbp --pmn 380 --weight 72
    python cli.py hrs --baseline-cr 1.0 --current-cr 2.3 --weight 70
    python cli.py tips --meld 16 --bili 2.2 --cr 1.1 --inr 1.3
    python cli.py batch --input sample.csv --output results.csv
    python cli.py interactive
"""

import argparse
import json
import os
import sys
from dataclasses import asdict

from cirrhosis_decompensation import (
    CirrhosisDecompensationEngine,
    AscitesDegree,
    EncephalopathyGrade,
    ChildPughClass,
    ACLFGrade,
    evaluate_meld_suite,
    calculate_child_pugh,
    evaluate_clif_aclf,
    evaluate_sbp_protocol,
    evaluate_hrs_aki_protocol,
    evaluate_tips_eligibility,
    process_batch_csv,
)


def format_cirrhosis_dossier_display(dossier):
    d = asdict(dossier)
    print("=" * 80)
    print(f"  CIRRHOSIS DECOMPENSATION & ACLF CLINICAL DOSSIER [{d['case_id']}]")
    print("=" * 80)
    print(f"  Patient ID:       {d['patient_id']}")
    print(f"  Timestamp (UTC):  {d['timestamp_utc']}")
    print("-" * 80)

    m = d['meld_suite']
    print(f"  MELD SUITE:")
    print(f"    - Original MELD (2002):  {m['original_meld']}")
    print(f"    - MELD-Na (UNOS 2016):   {m['meld_na']} (3-Month Waitlist Mortality: {m['three_month_mortality_pct']}%)")
    if m.get('meld_3_0'):
        print(f"    - MELD 3.0 (OPTN 2023):  {m['meld_3_0']}")
    print(f"    - Allocation Tier:       [{m['allocation_tier']}]")
    print("-" * 80)

    ctp = d['child_pugh']
    print(f"  CHILD-TURCOTTE-PUGH (CTP):")
    print(f"    - Class & Score:         Class {ctp['ctp_class']} ({ctp['total_points']} points)")
    print(f"    - Survival Estimates:    1-Year: {ctp['one_year_survival_pct']}% | 2-Year: {ctp['two_year_survival_pct']}%")
    print(f"    - Perioperative Risk:    {ctp['perioperative_mortality_pct']}% Mortality")
    print(f"    - Interpretation:        {ctp['clinical_interpretation']}")
    print("-" * 80)

    aclf = d['aclf_status']
    print(f"  EASL-CLIF ACLF STATUS:")
    print(f"    - Staging / Grade:       {aclf['aclf_grade_label']}")
    print(f"    - 28-Day Mortality:      {aclf['twenty_eight_day_mortality_pct']}%")
    print(f"    - ICU Indicated:         {'YES' if aclf['icu_admission_indicated'] else 'No'}")
    of = aclf['organ_failures']
    print(f"    - Organ Failures ({of['total_failures_count']}/6): "
          f"Liver={of['liver_failure']}, Kidney={of['kidney_failure']}, Brain={of['brain_failure']}, "
          f"Coagulation={of['coagulation_failure']}, Circulation={of['circulatory_failure']}, Respiration={of['respiratory_failure']}")
    print(f"    - Urgency Guidance:      {aclf['clinical_management_urgency']}")
    print("-" * 80)

    if d.get('sbp_protocol'):
        sbp = d['sbp_protocol']
        print(f"  SPONTANEOUS BACTERIAL PERITONITIS (SBP):")
        print(f"    - SBP Confirmed:         {'POSITIVE (PMN >= 250/mm³)' if sbp['is_sbp_confirmed'] else 'Negative'}")
        print(f"    - Ascitic PMN:           {sbp['ascitic_pmn_count']:.0f} / mm³")
        print(f"    - Antibiotic Regimen:    {sbp['antibiotic_regimen']}")
        if sbp['is_sbp_confirmed']:
            sch = sbp['albumin_dosing_schedule']
            print(f"    - Sort Albumin Infusion: Day 1 (within 6h): {sch['day_1_grams']} g | Day 3: {sch['day_3_grams']} g (Total: {sch['total_albumin_grams']} g)")
        print("-" * 80)

    if d.get('hrs_aki_protocol'):
        hrs = d['hrs_aki_protocol']
        print(f"  HEPATORENAL SYNDROME (HRS-AKI):")
        print(f"    - HRS-AKI Suspected:     {'YES (ICA Criteria Met)' if hrs['is_hrs_aki_suspected'] else 'No'}")
        print(f"    - KDIGO AKI Stage:       Stage {hrs['kdigo_aki_stage']}")
        print(f"    - First-Line Therapy:    {hrs['first_line_pharmacotherapy']}")
        print(f"    - Albumin Protocol:      {hrs['albumin_infusion_plan']}")
        print("-" * 80)

    tips = d['tips_eligibility']
    print(f"  TIPS ELIGIBILITY & SAFETY AUDIT:")
    print(f"    - Candidacy Status:      {'CANDIDATE' if tips['is_candidate'] else 'NON-CANDIDATE'}")
    print(f"    - Safety Risk Tier:      [{tips['risk_level']}]")
    if tips['absolute_contraindications']:
        print("    - Absolute Contraindications:")
        for c in tips['absolute_contraindications']:
            print(f"        ! {c}")
    if tips['relative_contraindications']:
        print("    - Relative Contraindications:")
        for c in tips['relative_contraindications']:
            print(f"        * {c}")
    print("-" * 80)

    if d['clinical_alerts']:
        print("  CRITICAL CLINICAL ALERTS:")
        for a in d['clinical_alerts']:
            print(f"    [!] {a}")
        print("-" * 80)

    print("  PRIORITIZED CLINICAL ACTION CHECKLIST:")
    for act in d['priority_action_checklist']:
        print(f"    {act}")
    print("=" * 80)


def cmd_evaluate(args):
    engine = CirrhosisDecompensationEngine()

    asc_map = {
        "none": AscitesDegree.NONE,
        "mild": AscitesDegree.MILD_CONTROLLED,
        "moderate": AscitesDegree.MODERATE_SEVERE_REFRACTORY,
        "severe": AscitesDegree.MODERATE_SEVERE_REFRACTORY,
        "refractory": AscitesDegree.MODERATE_SEVERE_REFRACTORY,
    }
    asc_val = asc_map.get(args.ascites.lower(), AscitesDegree.NONE)
    he_val = EncephalopathyGrade(args.he)

    dossier = engine.evaluate_patient_case(
        serum_creatinine_mg_dl=args.cr,
        total_bilirubin_mg_dl=args.bili,
        inr=args.inr,
        serum_sodium_mmol_l=args.na,
        serum_albumin_g_dl=args.alb,
        patient_weight_kg=args.weight,
        is_female=args.female,
        on_dialysis=args.dialysis,
        ascites=asc_val,
        encephalopathy=he_val,
        ascitic_pmn_count=args.pmn,
        baseline_creatinine_mg_dl=args.baseline_cr,
        requires_vasopressors=args.vasopressors,
        pao2_fio2_ratio=args.pf_ratio,
        has_severe_pulm_htn=args.pulm_htn,
        has_severe_heart_failure=args.heart_failure,
        case_id=args.case_id or "CASE-CIRR-001",
        patient_id=args.patient_id or "PT-HEP-001",
    )

    if args.json:
        print(json.dumps(asdict(dossier), indent=2, default=str))
    else:
        format_cirrhosis_dossier_display(dossier)
    return 0


def cmd_meld(args):
    result = evaluate_meld_suite(
        serum_creatinine_mg_dl=args.cr,
        total_bilirubin_mg_dl=args.bili,
        inr=args.inr,
        serum_sodium_mmol_l=args.na,
        serum_albumin_g_dl=args.alb,
        is_female=args.female,
        on_dialysis=args.dialysis,
    )

    if args.json:
        print(json.dumps(asdict(result), indent=2, default=str))
    else:
        print("=" * 60)
        print("  MELD SCORE SUITE (UNOS / OPTN GUIDELINES)")
        print("=" * 60)
        print(f"  Original MELD (2002):   {result.original_meld}")
        print(f"  MELD-Na (UNOS 2016):    {result.meld_na}")
        if result.meld_3_0:
            print(f"  MELD 3.0 (OPTN 2023):   {result.meld_3_0}")
        print(f"  3-Month Mortality:      {result.three_month_mortality_pct}%")
        print(f"  Allocation Tier:        {result.allocation_tier}")
        print(f"  Transplant Evaluation:  {'INDICATED (MELD-Na >= 15)' if result.details['transplant_evaluation_indicated'] else 'Monitor'}")
        print("=" * 60)
    return 0


def cmd_child_pugh(args):
    asc_map = {
        "none": AscitesDegree.NONE,
        "mild": AscitesDegree.MILD_CONTROLLED,
        "moderate": AscitesDegree.MODERATE_SEVERE_REFRACTORY,
        "severe": AscitesDegree.MODERATE_SEVERE_REFRACTORY,
    }
    asc_val = asc_map.get(args.ascites.lower(), AscitesDegree.NONE)
    he_val = EncephalopathyGrade(args.he)

    result = calculate_child_pugh(
        total_bilirubin_mg_dl=args.bili,
        serum_albumin_g_dl=args.alb,
        inr=args.inr,
        ascites=asc_val,
        encephalopathy=he_val,
        is_cholestatic_disease=args.cholestatic,
    )

    if args.json:
        print(json.dumps(asdict(result), indent=2, default=str))
    else:
        print("=" * 60)
        print("  CHILD-TURCOTTE-PUGH (CTP) SCORE & CLASSIFICATION")
        print("=" * 60)
        print(f"  Total Score:            {result.total_points} Points")
        print(f"  CTP Class:              Class {result.ctp_class.value}")
        print(f"  1-Year Survival:        {result.one_year_survival_pct}%")
        print(f"  2-Year Survival:        {result.two_year_survival_pct}%")
        print(f"  Perioperative Mortality:{result.perioperative_mortality_pct}%")
        print(f"\n  Clinical Summary:")
        print(f"  {result.clinical_interpretation}")
        print("=" * 60)
    return 0


def cmd_aclf(args):
    he_val = EncephalopathyGrade(args.he)
    result = evaluate_clif_aclf(
        total_bilirubin_mg_dl=args.bili,
        serum_creatinine_mg_dl=args.cr,
        inr=args.inr,
        encephalopathy_grade=he_val,
        requires_vasopressors=args.vasopressors,
        on_rrt_or_dialysis=args.dialysis,
        pao2_fio2_ratio=args.pf_ratio,
    )

    if args.json:
        print(json.dumps(asdict(result), indent=2, default=str))
    else:
        print("=" * 60)
        print("  EASL-CLIF ACLF STAGING & ORGAN FAILURE AUDIT")
        print("=" * 60)
        print(f"  ACLF Staging:           {result.aclf_grade_label}")
        print(f"  28-Day Mortality:       {result.twenty_eight_day_mortality_pct}%")
        print(f"  ICU Admission:          {'MANDATORY' if result.icu_admission_indicated else 'Ward / Step-down'}")
        print(f"\n  Clinical Management:")
        print(f"  {result.clinical_management_urgency}")
        print("=" * 60)
    return 0


def cmd_sbp(args):
    result = evaluate_sbp_protocol(
        ascitic_pmn_count_per_mm3=args.pmn,
        patient_weight_kg=args.weight,
    )

    if args.json:
        print(json.dumps(asdict(result), indent=2, default=str))
    else:
        print("=" * 60)
        print("  SPONTANEOUS BACTERIAL PERITONITIS (SBP) PROTOCOL")
        print("=" * 60)
        print(f"  Ascitic PMN Count:      {result.ascitic_pmn_count:.0f} / mm³")
        print(f"  Diagnosis:              {'POSITIVE SBP (PMN >= 250/mm³)' if result.is_sbp_confirmed else 'Negative (< 250/mm³)'}")
        print(f"  Antibiotic Regimen:     {result.antibiotic_regimen}")
        if result.is_sbp_confirmed:
            sch = result.albumin_dosing_schedule
            print(f"  IV Albumin Schedule:")
            print(f"    - Day 1 (1.5 g/kg):   {sch['day_1_grams']} g (within 6h)")
            print(f"    - Day 3 (1.0 g/kg):   {sch['day_3_grams']} g")
            print(f"    - Total Dose:         {sch['total_albumin_grams']} g")
        print("\n  Recommendations:")
        for r in result.recommendations:
            print(f"    * {r}")
        print("=" * 60)
    return 0


def cmd_hrs(args):
    result = evaluate_hrs_aki_protocol(
        baseline_creatinine_mg_dl=args.baseline_cr,
        current_creatinine_mg_dl=args.current_cr,
        patient_weight_kg=args.weight,
        has_ascites=not args.no_ascites,
    )

    if args.json:
        print(json.dumps(asdict(result), indent=2, default=str))
    else:
        print("=" * 60)
        print("  HEPATORENAL SYNDROME (HRS-AKI) ICA EVALUATION")
        print("=" * 60)
        print(f"  HRS-AKI Diagnosis:      {'CONFIRMED' if result.is_hrs_aki_suspected else 'Not Met'}")
        print(f"  KDIGO AKI Staging:      Stage {result.kdigo_aki_stage}")
        print(f"  Pharmacotherapy:        {result.first_line_pharmacotherapy}")
        print(f"  Albumin Protocol:       {result.albumin_infusion_plan}")
        print("\n  Clinical Guidance:")
        for g in result.management_guidance:
            print(f"    * {g}")
        print("=" * 60)
    return 0


def cmd_tips(args):
    result = evaluate_tips_eligibility(
        meld_score=args.meld,
        total_bilirubin_mg_dl=args.bili,
        serum_creatinine_mg_dl=args.cr,
        inr=args.inr,
        has_severe_pulmonary_hypertension=args.pulm_htn,
        has_congestive_heart_failure=args.heart_failure,
        has_recurrent_severe_encephalopathy=args.encephalopathy,
    )

    if args.json:
        print(json.dumps(asdict(result), indent=2, default=str))
    else:
        print("=" * 60)
        print("  TIPS ELIGIBILITY & PRE-PROCEDURAL SAFETY AUDIT")
        print("=" * 60)
        print(f"  Candidacy:              {'ELIGIBLE' if result.is_candidate else 'CONTRAINDICATED'}")
        print(f"  Risk Classification:    [{result.risk_level}]")
        if result.absolute_contraindications:
            print("\n  Absolute Contraindications:")
            for c in result.absolute_contraindications:
                print(f"    ! {c}")
        if result.relative_contraindications:
            print("\n  Relative Contraindications:")
            for c in result.relative_contraindications:
                print(f"    * {c}")
        print("\n  Recommendations:")
        for r in result.recommendations:
            print(f"    - {r}")
        print("=" * 60)
    return 0


def cmd_batch(args):
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found.", file=sys.stderr)
        return 1
    count = process_batch_csv(args.input, args.output)
    print(f"Successfully processed {count} cirrhosis records from '{args.input}' -> '{args.output}'.")
    return 0


def cmd_interactive(args):
    print("=" * 70)
    print("  CIRRHOSIS DECOMPENSATION & ACLF INTERACTIVE CLINICAL WIZARD")
    print("=" * 70)
    cr = float(input("Serum Creatinine (mg/dL, default 1.2): ") or "1.2")
    bili = float(input("Total Bilirubin (mg/dL, default 2.5): ") or "2.5")
    inr = float(input("INR (default 1.5): ") or "1.5")
    na = float(input("Serum Sodium (mmol/L, default 133): ") or "133")
    alb = float(input("Serum Albumin (g/dL, default 2.8): ") or "2.8")
    wt = float(input("Patient Weight (kg, default 70): ") or "70")
    female = input("Is patient female? (y/n, default n): ").lower().startswith("y")
    dialysis = input("On renal replacement / hemodialysis? (y/n, default n): ").lower().startswith("y")
    
    print("\nAscites Status: [1] None  [2] Mild/Controlled  [3] Moderate/Severe/Refractory")
    a_ch = input("Choice (default 2): ") or "2"
    asc_map = {"1": AscitesDegree.NONE, "2": AscitesDegree.MILD_CONTROLLED, "3": AscitesDegree.MODERATE_SEVERE_REFRACTORY}
    asc = asc_map.get(a_ch, AscitesDegree.MILD_CONTROLLED)

    he_val = int(input("\nHepatic Encephalopathy Grade (0 to 4, default 1): ") or "1")
    pmn_str = input("Ascitic PMN count / mm³ (optional, e.g. 300): ")
    pmn = float(pmn_str) if pmn_str else None

    engine = CirrhosisDecompensationEngine()
    dossier = engine.evaluate_patient_case(
        serum_creatinine_mg_dl=cr,
        total_bilirubin_mg_dl=bili,
        inr=inr,
        serum_sodium_mmol_l=na,
        serum_albumin_g_dl=alb,
        patient_weight_kg=wt,
        is_female=female,
        on_dialysis=dialysis,
        ascites=asc,
        encephalopathy=EncephalopathyGrade(he_val) if he_val in range(5) else EncephalopathyGrade.GRADE_0,
        ascitic_pmn_count=pmn,
    )

    print("\n")
    format_cirrhosis_dossier_display(dossier)
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="cirrhosis-decompensation-agent",
        description="Cirrhosis Decompensation, MELD Suite & ACLF Clinical Decision Support Engine",
    )
    parser.add_argument("--json", action="store_true", help="Output results in JSON format")

    subparsers = parser.add_subparsers(dest="command")

    # Evaluate
    p_eval = subparsers.add_parser("evaluate", help="Comprehensive multi-system cirrhosis audit")
    p_eval.add_argument("--cr", type=float, required=True, help="Serum Creatinine (mg/dL)")
    p_eval.add_argument("--bili", type=float, required=True, help="Total Bilirubin (mg/dL)")
    p_eval.add_argument("--inr", type=float, required=True, help="INR")
    p_eval.add_argument("--na", type=float, default=135.0, help="Serum Sodium (mmol/L)")
    p_eval.add_argument("--alb", type=float, default=3.0, help="Serum Albumin (g/dL)")
    p_eval.add_argument("--weight", type=float, default=70.0, help="Weight (kg)")
    p_eval.add_argument("--female", action="store_true", help="Female patient")
    p_eval.add_argument("--dialysis", action="store_true", help="On dialysis / RRT")
    p_eval.add_argument("--ascites", default="none", choices=["none", "mild", "moderate", "severe", "refractory"])
    p_eval.add_argument("--he", type=int, default=0, choices=[0, 1, 2, 3, 4], help="West Haven HE grade")
    p_eval.add_argument("--pmn", type=float, help="Ascitic PMN count / mm³")
    p_eval.add_argument("--baseline-cr", type=float, help="Baseline Creatinine (mg/dL)")
    p_eval.add_argument("--vasopressors", action="store_true", help="Requiring vasopressors")
    p_eval.add_argument("--pf-ratio", type=float, help="PaO2/FiO2 ratio")
    p_eval.add_argument("--pulm-htn", action="store_true", help="Severe pulmonary HTN")
    p_eval.add_argument("--heart-failure", action="store_true", help="Severe heart failure")
    p_eval.add_argument("--case-id", help="Case ID")
    p_eval.add_argument("--patient-id", help="Patient ID")

    # MELD
    p_meld = subparsers.add_parser("meld", help="Calculate MELD, MELD-Na, and MELD 3.0")
    p_meld.add_argument("--cr", type=float, required=True, help="Serum Creatinine (mg/dL)")
    p_meld.add_argument("--bili", type=float, required=True, help="Total Bilirubin (mg/dL)")
    p_meld.add_argument("--inr", type=float, required=True, help="INR")
    p_meld.add_argument("--na", type=float, default=135.0, help="Serum Sodium (mmol/L)")
    p_meld.add_argument("--alb", type=float, help="Serum Albumin (g/dL, for MELD 3.0)")
    p_meld.add_argument("--female", action="store_true", help="Female gender (for MELD 3.0)")
    p_meld.add_argument("--dialysis", action="store_true", help="On dialysis >= 2x in past 7 days")

    # Child-Pugh
    p_ctp = subparsers.add_parser("child-pugh", help="Calculate Child-Turcotte-Pugh score & class")
    p_ctp.add_argument("--bili", type=float, required=True, help="Total Bilirubin (mg/dL)")
    p_ctp.add_argument("--alb", type=float, required=True, help="Serum Albumin (g/dL)")
    p_ctp.add_argument("--inr", type=float, required=True, help="INR")
    p_ctp.add_argument("--ascites", default="none", choices=["none", "mild", "moderate", "severe"])
    p_ctp.add_argument("--he", type=int, default=0, choices=[0, 1, 2, 3, 4], help="West Haven HE grade")
    p_ctp.add_argument("--cholestatic", action="store_true", help="Cholestatic liver disease (PBC/PSC)")

    # ACLF
    p_aclf = subparsers.add_parser("aclf", help="Evaluate EASL-CLIF ACLF staging")
    p_aclf.add_argument("--cr", type=float, required=True, help="Serum Creatinine (mg/dL)")
    p_aclf.add_argument("--bili", type=float, required=True, help="Total Bilirubin (mg/dL)")
    p_aclf.add_argument("--inr", type=float, required=True, help="INR")
    p_aclf.add_argument("--he", type=int, default=0, choices=[0, 1, 2, 3, 4], help="West Haven HE grade")
    p_aclf.add_argument("--vasopressors", action="store_true", help="Vasopressors required")
    p_aclf.add_argument("--dialysis", action="store_true", help="On RRT / dialysis")
    p_aclf.add_argument("--pf-ratio", type=float, help="PaO2 / FiO2 ratio")

    # SBP
    p_sbp = subparsers.add_parser("sbp", help="Evaluate SBP criteria & Sort Albumin protocol")
    p_sbp.add_argument("--pmn", type=float, required=True, help="Ascitic PMN count / mm³")
    p_sbp.add_argument("--weight", type=float, required=True, help="Patient weight in kg")

    # HRS
    p_hrs = subparsers.add_parser("hrs", help="Evaluate HRS-AKI criteria & Terlipressin regimen")
    p_hrs.add_argument("--baseline-cr", type=float, required=True, help="Baseline Creatinine (mg/dL)")
    p_hrs.add_argument("--current-cr", type=float, required=True, help="Current Creatinine (mg/dL)")
    p_hrs.add_argument("--weight", type=float, required=True, help="Weight in kg")
    p_hrs.add_argument("--no-ascites", action="store_true", help="Set if ascites absent")

    # TIPS
    p_tips = subparsers.add_parser("tips", help="Evaluate TIPS eligibility and contraindications")
    p_tips.add_argument("--meld", type=int, required=True, help="MELD-Na score")
    p_tips.add_argument("--bili", type=float, required=True, help="Total Bilirubin (mg/dL)")
    p_tips.add_argument("--cr", type=float, required=True, help="Serum Creatinine (mg/dL)")
    p_tips.add_argument("--inr", type=float, required=True, help="INR")
    p_tips.add_argument("--pulm-htn", action="store_true", help="Severe pulmonary hypertension")
    p_tips.add_argument("--heart-failure", action="store_true", help="Severe heart failure")
    p_tips.add_argument("--encephalopathy", action="store_true", help="Recurrent severe encephalopathy")

    # Batch
    p_batch = subparsers.add_parser("batch", help="Batch process cirrhosis cases from CSV")
    p_batch.add_argument("-i", "--input", required=True, help="Input CSV file")
    p_batch.add_argument("-o", "--output", default="cirrhosis_results.csv", help="Output CSV file")

    # Interactive
    p_inter = subparsers.add_parser("interactive", help="Interactive clinical wizard")

    args = parser.parse_args(argv)

    if args.command == "evaluate":
        return cmd_evaluate(args)
    elif args.command == "meld":
        return cmd_meld(args)
    elif args.command == "child-pugh":
        return cmd_child_pugh(args)
    elif args.command == "aclf":
        return cmd_aclf(args)
    elif args.command == "sbp":
        return cmd_sbp(args)
    elif args.command == "hrs":
        return cmd_hrs(args)
    elif args.command == "tips":
        return cmd_tips(args)
    elif args.command == "batch":
        return cmd_batch(args)
    elif args.command == "interactive":
        return cmd_interactive(args)
    else:
        parser.print_help()
        return 0


if __name__ == "__main__":
    sys.exit(main())
