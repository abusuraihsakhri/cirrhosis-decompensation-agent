#!/usr/bin/env python3
"""
Unit Test Suite for Cirrhosis Decompensation & ACLF Clinical Decision Support Engine.
"""

import os
import tempfile
import unittest
from cirrhosis_decompensation import (
    CirrhosisDecompensationEngine,
    AscitesDegree,
    EncephalopathyGrade,
    ChildPughClass,
    ACLFGrade,
    calculate_original_meld,
    calculate_meld_na,
    calculate_meld_3_0,
    calculate_child_pugh,
    evaluate_clif_aclf,
    evaluate_sbp_protocol,
    evaluate_hrs_aki_protocol,
    evaluate_tips_eligibility,
    evaluate_meld_suite,
    process_batch_csv,
)


class TestMELDCalculations(unittest.TestCase):
    """Test original MELD, MELD-Na, and MELD 3.0 calculation accuracy."""

    def test_original_meld_standard(self):
        # Cr=1.5, Bili=2.5, INR=1.3 -> 9.57*ln(1.5) + 3.78*ln(2.5) + 11.20*ln(1.3) + 6.43
        # = 9.57*0.40546 + 3.78*0.91629 + 11.20*0.26236 + 6.43
        # = 3.880 + 3.464 + 2.938 + 6.43 = 16.712 -> 17
        score = calculate_original_meld(serum_creatinine_mg_dl=1.5, total_bilirubin_mg_dl=2.5, inr=1.3)
        self.assertEqual(score, 17)

    def test_original_meld_lower_bounding(self):
        # All values below 1.0 bounded up to 1.0 -> ln(1) = 0 -> score = 6.43 -> 6
        score = calculate_original_meld(serum_creatinine_mg_dl=0.6, total_bilirubin_mg_dl=0.7, inr=0.9)
        self.assertEqual(score, 6)

    def test_original_meld_creatinine_capping(self):
        # Cr=6.0 (capped at 4.0), Bili=1.0, INR=1.0 -> 9.57*ln(4.0) + 6.43 = 9.57*1.38629 + 6.43 = 13.267 + 6.43 = 19.697 -> 20
        score = calculate_original_meld(serum_creatinine_mg_dl=6.0, total_bilirubin_mg_dl=1.0, inr=1.0)
        self.assertEqual(score, 20)

    def test_original_meld_dialysis_override(self):
        # On dialysis -> Cr automatically set to 4.0
        score = calculate_original_meld(serum_creatinine_mg_dl=1.2, total_bilirubin_mg_dl=1.0, inr=1.0, on_dialysis=True)
        self.assertEqual(score, 20)

    def test_meld_na_when_meld_le_11(self):
        # If MELD <= 11, MELD-Na equals original MELD regardless of sodium
        meld_na = calculate_meld_na(original_meld=9, serum_sodium_mmol_l=125.0)
        self.assertEqual(meld_na, 9)

    def test_meld_na_hyponatremia_adjustment(self):
        # MELD=20, Na=128 (137-128 = 9)
        # MELD-Na = 20 + 1.32*9 - [0.033 * 20 * 9] = 20 + 11.88 - 5.94 = 25.94 -> 26
        meld_na = calculate_meld_na(original_meld=20, serum_sodium_mmol_l=128.0)
        self.assertEqual(meld_na, 26)

    def test_meld_na_sodium_bounds(self):
        # Sodium < 125 bounded to 125
        meld_na_120 = calculate_meld_na(original_meld=20, serum_sodium_mmol_l=120.0)
        meld_na_125 = calculate_meld_na(original_meld=20, serum_sodium_mmol_l=125.0)
        self.assertEqual(meld_na_120, meld_na_125)

    def test_meld_3_0_female_bonus_and_albumin(self):
        # Evaluate MELD 3.0 calculation for female vs male
        score_f = calculate_meld_3_0(
            serum_creatinine_mg_dl=1.5,
            total_bilirubin_mg_dl=2.5,
            inr=1.3,
            serum_sodium_mmol_l=133.0,
            serum_albumin_g_dl=2.8,
            is_female=True,
        )
        score_m = calculate_meld_3_0(
            serum_creatinine_mg_dl=1.5,
            total_bilirubin_mg_dl=2.5,
            inr=1.3,
            serum_sodium_mmol_l=133.0,
            serum_albumin_g_dl=2.8,
            is_female=False,
        )
        self.assertGreaterEqual(score_f, score_m)

    def test_meld_invalid_inputs(self):
        with self.assertRaises(ValueError):
            calculate_original_meld(-1.0, 2.0, 1.0)
        with self.assertRaises(ValueError):
            calculate_meld_na(15, -130.0)


class TestChildPughScoring(unittest.TestCase):
    """Test Child-Turcotte-Pugh score and classification."""

    def test_class_a_compensated(self):
        # Bili 1.2 (1), Alb 3.8 (1), INR 1.1 (1), Ascites None (1), HE None (1) -> 5 pts (Class A)
        res = calculate_child_pugh(
            total_bilirubin_mg_dl=1.2,
            serum_albumin_g_dl=3.8,
            inr=1.1,
            ascites=AscitesDegree.NONE,
            encephalopathy=EncephalopathyGrade.GRADE_0,
        )
        self.assertEqual(res.total_points, 5)
        self.assertEqual(res.ctp_class, ChildPughClass.CLASS_A)
        self.assertEqual(res.one_year_survival_pct, 100.0)

    def test_class_b_moderate(self):
        # Bili 2.5 (2), Alb 3.0 (2), INR 1.9 (2), Ascites Mild (2), HE None (1) -> 9 pts (Class B)
        res = calculate_child_pugh(
            total_bilirubin_mg_dl=2.5,
            serum_albumin_g_dl=3.0,
            inr=1.9,
            ascites=AscitesDegree.MILD_CONTROLLED,
            encephalopathy=EncephalopathyGrade.GRADE_0,
        )
        self.assertEqual(res.total_points, 9)
        self.assertEqual(res.ctp_class, ChildPughClass.CLASS_B)
        self.assertEqual(res.one_year_survival_pct, 80.0)

    def test_class_c_decompensated(self):
        # Bili 4.5 (3), Alb 2.2 (3), INR 2.5 (3), Ascites Moderate (3), HE Grade 3 (3) -> 15 pts (Class C)
        res = calculate_child_pugh(
            total_bilirubin_mg_dl=4.5,
            serum_albumin_g_dl=2.2,
            inr=2.5,
            ascites=AscitesDegree.MODERATE_SEVERE_REFRACTORY,
            encephalopathy=EncephalopathyGrade.GRADE_3,
        )
        self.assertEqual(res.total_points, 15)
        self.assertEqual(res.ctp_class, ChildPughClass.CLASS_C)
        self.assertEqual(res.one_year_survival_pct, 45.0)

    def test_cholestatic_disease_thresholds(self):
        # In cholestatic disease (PBC/PSC), Bili < 4.0 is 1 pt (instead of < 2.0)
        res_normal = calculate_child_pugh(
            total_bilirubin_mg_dl=3.5,
            serum_albumin_g_dl=4.0,
            inr=1.0,
            ascites=AscitesDegree.NONE,
            encephalopathy=EncephalopathyGrade.GRADE_0,
            is_cholestatic_disease=False,
        )
        res_cholestatic = calculate_child_pugh(
            total_bilirubin_mg_dl=3.5,
            serum_albumin_g_dl=4.0,
            inr=1.0,
            ascites=AscitesDegree.NONE,
            encephalopathy=EncephalopathyGrade.GRADE_0,
            is_cholestatic_disease=True,
        )
        self.assertEqual(res_normal.point_breakdown["bilirubin_points"], 3)
        self.assertEqual(res_cholestatic.point_breakdown["bilirubin_points"], 1)


class TestEASLCLIFAssessment(unittest.TestCase):
    """Test ACLF organ failures and staging."""

    def test_no_aclf(self):
        res = evaluate_clif_aclf(
            total_bilirubin_mg_dl=3.0,
            serum_creatinine_mg_dl=1.1,
            inr=1.3,
            encephalopathy_grade=EncephalopathyGrade.GRADE_0,
        )
        self.assertEqual(res.aclf_grade, ACLFGrade.NO_ACLF)
        self.assertFalse(res.icu_admission_indicated)

    def test_aclf_grade_1_single_kidney_failure(self):
        res = evaluate_clif_aclf(
            total_bilirubin_mg_dl=3.0,
            serum_creatinine_mg_dl=2.4,  # >= 2.0 Kidney Failure
            inr=1.4,
            encephalopathy_grade=EncephalopathyGrade.GRADE_0,
        )
        self.assertEqual(res.aclf_grade, ACLFGrade.GRADE_1)
        self.assertTrue(res.organ_failures.kidney_failure)
        self.assertEqual(res.organ_failures.total_failures_count, 1)

    def test_aclf_grade_2_two_organ_failures(self):
        res = evaluate_clif_aclf(
            total_bilirubin_mg_dl=14.0,  # Liver failure (>=12)
            serum_creatinine_mg_dl=2.5,  # Kidney failure (>=2.0)
            inr=1.5,
            encephalopathy_grade=EncephalopathyGrade.GRADE_1,
        )
        self.assertEqual(res.aclf_grade, ACLFGrade.GRADE_2)
        self.assertTrue(res.icu_admission_indicated)
        self.assertEqual(res.organ_failures.total_failures_count, 2)

    def test_aclf_grade_3_three_or_more_failures(self):
        res = evaluate_clif_aclf(
            total_bilirubin_mg_dl=15.0,  # Liver failure
            serum_creatinine_mg_dl=2.8,  # Kidney failure
            inr=2.8,                     # Coagulation failure (>=2.5)
            encephalopathy_grade=EncephalopathyGrade.GRADE_3,  # Brain failure
            requires_vasopressors=True,  # Circulatory failure
        )
        self.assertEqual(res.aclf_grade, ACLFGrade.GRADE_3)
        self.assertEqual(res.organ_failures.total_failures_count, 5)
        self.assertEqual(res.twenty_eight_day_mortality_pct, 78.6)


class TestAcuteDecompensationProtocols(unittest.TestCase):
    """Test SBP, HRS-AKI, and TIPS protocol engines."""

    def test_sbp_positive_sort_albumin_dosing(self):
        # PMN 350 (>=250) in 70kg patient -> SBP positive
        # Day 1: 1.5 * 70 = 105.0 g, Day 3: 1.0 * 70 = 70.0 g (Total = 175.0 g)
        res = evaluate_sbp_protocol(ascitic_pmn_count_per_mm3=350.0, patient_weight_kg=70.0)
        self.assertTrue(res.is_sbp_confirmed)
        self.assertEqual(res.albumin_dosing_schedule["day_1_grams"], 105.0)
        self.assertEqual(res.albumin_dosing_schedule["day_3_grams"], 70.0)
        self.assertEqual(res.albumin_dosing_schedule["total_albumin_grams"], 175.0)
        self.assertIn("Ceftriaxone", res.antibiotic_regimen)

    def test_sbp_negative(self):
        res = evaluate_sbp_protocol(ascitic_pmn_count_per_mm3=120.0, patient_weight_kg=70.0)
        self.assertFalse(res.is_sbp_confirmed)
        self.assertEqual(res.albumin_dosing_schedule["day_1_grams"], 0)

    def test_hrs_aki_confirmed(self):
        res = evaluate_hrs_aki_protocol(
            baseline_creatinine_mg_dl=0.9,
            current_creatinine_mg_dl=2.2,  # >2x baseline -> KDIGO Stage 2
            patient_weight_kg=70.0,
            has_ascites=True,
            no_response_to_48h_albumin_expansion=True,
        )
        self.assertTrue(res.is_hrs_aki_suspected)
        self.assertEqual(res.kdigo_aki_stage, 2)
        self.assertIn("Terlipressin", res.first_line_pharmacotherapy)

    def test_tips_contraindicated_pulmonary_hypertension(self):
        res = evaluate_tips_eligibility(
            meld_score=14,
            total_bilirubin_mg_dl=1.8,
            serum_creatinine_mg_dl=1.0,
            inr=1.2,
            has_severe_pulmonary_hypertension=True,
        )
        self.assertFalse(res.is_candidate)
        self.assertEqual(res.risk_level, "CONTRAINDICATED")
        self.assertTrue(len(res.absolute_contraindications) >= 1)

    def test_tips_favorable_candidate(self):
        res = evaluate_tips_eligibility(
            meld_score=12,
            total_bilirubin_mg_dl=1.5,
            serum_creatinine_mg_dl=0.9,
            inr=1.2,
        )
        self.assertTrue(res.is_candidate)
        self.assertEqual(res.risk_level, "FAVORABLE CANDIDATE")


class TestMasterEngineAndBatch(unittest.TestCase):
    """Test full engine orchestration and CSV processing."""

    def setUp(self):
        self.engine = CirrhosisDecompensationEngine()

    def test_full_patient_dossier_generation(self):
        dossier = self.engine.evaluate_patient_case(
            serum_creatinine_mg_dl=1.8,
            total_bilirubin_mg_dl=3.2,
            inr=1.7,
            serum_sodium_mmol_l=131.0,
            serum_albumin_g_dl=2.6,
            patient_weight_kg=72.0,
            is_female=True,
            ascites=AscitesDegree.MODERATE_SEVERE_REFRACTORY,
            encephalopathy=EncephalopathyGrade.GRADE_2,
            ascitic_pmn_count=400.0,
            baseline_creatinine_mg_dl=1.0,
        )
        self.assertGreater(dossier.meld_suite.meld_na, 15)
        self.assertEqual(dossier.child_pugh.ctp_class, ChildPughClass.CLASS_C)
        self.assertTrue(dossier.sbp_protocol.is_sbp_confirmed)
        self.assertTrue(len(dossier.clinical_alerts) >= 2)
        self.assertTrue(len(dossier.priority_action_checklist) >= 3)

    def test_meld_3_0_interaction_terms(self):
        # High bilirubin and severe hyponatremia with interaction term
        score = calculate_meld_3_0(
            serum_creatinine_mg_dl=2.5,
            total_bilirubin_mg_dl=10.0,
            inr=2.2,
            serum_sodium_mmol_l=126.0,
            serum_albumin_g_dl=2.2,
            is_female=True,
            on_dialysis=False,
        )
        self.assertGreaterEqual(score, 30)
        self.assertLessEqual(score, 40)

    def test_meld_suite_status_1a_tier(self):
        res = evaluate_meld_suite(
            serum_creatinine_mg_dl=3.8,
            total_bilirubin_mg_dl=18.0,
            inr=3.2,
            serum_sodium_mmol_l=125.0,
            serum_albumin_g_dl=2.0,
            is_female=False,
        )
        self.assertGreaterEqual(res.meld_na, 35)
        self.assertIn("STATUS 1A", res.allocation_tier)
        self.assertTrue(res.details["transplant_evaluation_indicated"])

    def test_aclf_respiratory_failure_pao2_fio2(self):
        # PaO2/FiO2 <= 200 is respiratory failure
        res = evaluate_clif_aclf(
            total_bilirubin_mg_dl=2.0,
            serum_creatinine_mg_dl=1.0,
            inr=1.2,
            encephalopathy_grade=EncephalopathyGrade.GRADE_0,
            pao2_fio2_ratio=180.0,
        )
        self.assertTrue(res.organ_failures.respiratory_failure)

    def test_aclf_respiratory_failure_spo2_fio2(self):
        # SpO2/FiO2 <= 214 is respiratory failure
        res = evaluate_clif_aclf(
            total_bilirubin_mg_dl=2.0,
            serum_creatinine_mg_dl=1.0,
            inr=1.2,
            encephalopathy_grade=EncephalopathyGrade.GRADE_0,
            spo2_fio2_ratio=205.0,
        )
        self.assertTrue(res.organ_failures.respiratory_failure)

    def test_hrs_aki_stage_3_high_creatinine(self):
        res = evaluate_hrs_aki_protocol(
            baseline_creatinine_mg_dl=1.1,
            current_creatinine_mg_dl=4.2,  # Cr >= 4.0 -> Stage 3
            patient_weight_kg=75.0,
            has_ascites=True,
        )
        self.assertEqual(res.kdigo_aki_stage, 3)
        self.assertTrue(res.is_hrs_aki_suspected)

    def test_tips_contraindicated_congestive_heart_failure(self):
        res = evaluate_tips_eligibility(
            meld_score=13,
            total_bilirubin_mg_dl=1.2,
            serum_creatinine_mg_dl=0.9,
            inr=1.1,
            has_congestive_heart_failure=True,
        )
        self.assertFalse(res.is_candidate)
        self.assertEqual(res.risk_level, "CONTRAINDICATED")
    def test_batch_csv_processing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            in_csv = os.path.join(tmpdir, "cirrhosis_input.csv")
            out_csv = os.path.join(tmpdir, "cirrhosis_output.csv")

            with open(in_csv, "w", encoding="utf-8") as f:
                f.write("case_id,patient_id,creatinine,bilirubin,inr,sodium,albumin,weight_kg,is_female,ascites,encephalopathy_grade\n")
                f.write("CASE-001,PT-101,1.1,1.5,1.2,138,3.6,75,false,none,0\n")
                f.write("CASE-002,PT-102,2.2,4.0,2.1,128,2.5,65,true,moderate,2\n")

            count = process_batch_csv(in_csv, out_csv)
            self.assertEqual(count, 2)
            self.assertTrue(os.path.exists(out_csv))

            with open(out_csv, "r", encoding="utf-8") as f:
                lines = f.readlines()
                self.assertEqual(len(lines), 3)


if __name__ == "__main__":
    unittest.main()
