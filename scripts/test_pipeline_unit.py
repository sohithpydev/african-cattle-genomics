#!/usr/bin/env python3
"""
Unit Smoke Tests & Data Integrity Checks for Livestock Genomics Pipeline
========================================================================
Author: Sohith Reddy Korlagunta
Target: Reproducibility & Integrity Verification
"""

import unittest
import importlib.util
import os
import numpy as np
import pandas as pd

def load_module(name, file_path):
    spec = importlib.util.spec_from_file_location(name, file_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

mod_pop = load_module("pop_struct", "scripts/02_population_structure.py")
classify_lineage = mod_pop.classify_lineage

mod_roh = load_module("roh_detect", "scripts/03_detect_roh.py")
detect_roh_individual = mod_roh.detect_roh_individual
classify_length = mod_roh.classify_length

class TestLivestockPipeline(unittest.TestCase):

    def test_metadata_integrity(self):
        meta = pd.read_csv("data/metadata/sample_metadata.tsv", sep="\t")
        self.assertEqual(len(meta), 1082, "Total metadata records must be 1082.")
        self.assertEqual(sum(meta["Study_Cohort"] == "INCLUDED"), 1063, "Included samples must be exactly 1063.")
        self.assertEqual(sum(meta["Study_Cohort"] == "EXCLUDED"), 19, "Excluded samples must be exactly 19.")
        self.assertEqual(meta["Breed"].nunique(), 40, "Cohort must span exactly 40 breeds.")

    def test_lineage_classification(self):
        self.assertEqual(classify_lineage("West African Taurine", "Ndama"), "African Taurine")
        self.assertEqual(classify_lineage("East African Zebu", "Ethiopia Boran"), "African Zebu")
        self.assertEqual(classify_lineage("East African Admixed", "Ankole"), "Sanga / Zenga (Admixed)")
        self.assertEqual(classify_lineage("Central African Admixed", "Kuri"), "Sanga / Zenga (Admixed)")

    def test_roh_length_classification(self):
        self.assertEqual(classify_length(1.5), "1-2 Mb")
        self.assertEqual(classify_length(2.5), "2-4 Mb")
        self.assertEqual(classify_length(5.0), "4-8 Mb")
        self.assertEqual(classify_length(12.0), "8-16 Mb")
        self.assertEqual(classify_length(22.0), ">=16 Mb")

    def test_roh_detector_synthetic(self):
        # Create synthetic chromosome of 100 SNPs spanning 2.0 Mb
        n_snps = 100
        chrom_arr = np.ones(n_snps, dtype=int)
        pos_arr = np.linspace(1000000, 3000000, n_snps, dtype=int) # 2.0 Mb length, ~20kb spacing

        # Test 1: Completely homozygous segment
        calls_hom = np.zeros(n_snps, dtype=int) # all 0/0
        rohs = detect_roh_individual(calls_hom, chrom_arr, pos_arr, "TEST_01", min_length_bp=1000000, min_snps=50)
        self.assertEqual(len(rohs), 1, "Should detect exactly 1 ROH segment.")
        self.assertEqual(rohs[0]["Length_Class"], "2-4 Mb")

        # Test 2: Highly heterozygous segment (should detect 0)
        calls_het = np.random.choice([0, 1], size=n_snps, p=[0.5, 0.5])
        rohs_het = detect_roh_individual(calls_het, chrom_arr, pos_arr, "TEST_02", min_length_bp=1000000, min_snps=50)
        self.assertEqual(len(rohs_het), 0, "Heterozygous region should yield 0 ROH.")

if __name__ == "__main__":
    unittest.main()
