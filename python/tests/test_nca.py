"""
Tests for OpenPKPD NCA (Non-Compartmental Analysis) module.

These tests verify the Python bindings to the Julia NCA implementation.
"""

import pytest
import openpkpd
from openpkpd.nca import (
    run_nca,
    NCAConfig,
    NCAResult,
    nca_cmax,
    nca_tmax,
    nca_cmin,
    nca_clast,
    nca_tlast,
    auc_0_t,
    estimate_lambda_z,
    nca_half_life,
    nca_mrt,
    nca_cl_f,
    nca_vz_f,
    nca_accumulation_index,
    nca_ptf,
    nca_swing,
    bioequivalence_90ci,
    tost_analysis,
    be_conclusion,
    geometric_mean_ratio,
    geometric_mean,
    within_subject_cv,
)


# ============================================================================
# Test Data
# ============================================================================

# Standard oral PK profile
TEST_TIMES = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0]
TEST_CONC = [0.0, 0.82, 1.44, 1.62, 1.28, 0.94, 0.68, 0.36, 0.08]
TEST_DOSE = 100.0

# Multiple dose profile
SS_TIMES = [0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0]
SS_CONC = [2.5, 4.2, 3.8, 3.0, 2.8, 2.6, 2.5, 2.5]
SS_TAU = 12.0


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture(scope="module")
def init():
    """Initialize Julia once for all tests."""
    openpkpd.init_julia()


# ============================================================================
# Primary Exposure Metrics Tests
# ============================================================================

class TestExposureMetrics:
    """Tests for primary exposure metrics."""

    def test_nca_cmax(self, init):
        """Test Cmax calculation."""
        cmax = nca_cmax(TEST_CONC)
        assert cmax == 1.62

    def test_nca_tmax(self, init):
        """Test Tmax calculation."""
        tmax = nca_tmax(TEST_TIMES, TEST_CONC)
        assert tmax == 2.0

    def test_nca_cmin(self, init):
        """Test Cmin calculation."""
        cmin = nca_cmin(SS_CONC)
        assert cmin == 2.5

    def test_nca_clast(self, init):
        """Test Clast calculation."""
        clast = nca_clast(TEST_TIMES, TEST_CONC)
        assert clast == 0.08

    def test_nca_tlast(self, init):
        """Test Tlast calculation."""
        tlast = nca_tlast(TEST_TIMES, TEST_CONC)
        assert tlast == 24.0

    def test_clast_with_lloq(self, init):
        """Test Clast with LLOQ."""
        clast = nca_clast(TEST_TIMES, TEST_CONC, lloq=0.1)
        assert clast == 0.36

    def test_tlast_with_lloq(self, init):
        """Test Tlast with LLOQ."""
        tlast = nca_tlast(TEST_TIMES, TEST_CONC, lloq=0.1)
        assert tlast == 12.0


# ============================================================================
# Lambda-z and AUC Tests
# ============================================================================

class TestLambdaZAndAUC:
    """Tests for lambda-z estimation and AUC calculations."""

    def test_estimate_lambda_z(self, init):
        """Test lambda-z estimation."""
        result = estimate_lambda_z(TEST_TIMES, TEST_CONC)

        assert result["lambda_z"] is not None
        assert result["lambda_z"] > 0.0
        assert result["t_half"] is not None
        assert result["t_half"] > 0.0
        assert result["r_squared"] is not None
        assert result["r_squared"] >= 0.9
        assert result["quality_flag"] in ["good", "warning", "insufficient"]

    def test_lambda_z_t_half_relationship(self, init):
        """Test that t_half = ln(2) / lambda_z."""
        result = estimate_lambda_z(TEST_TIMES, TEST_CONC)

        if result["lambda_z"] is not None and result["t_half"] is not None:
            import math
            expected_t_half = math.log(2) / result["lambda_z"]
            assert abs(result["t_half"] - expected_t_half) < 1e-10

    def test_nca_half_life(self, init):
        """Test half-life calculation from lambda_z."""
        import math
        lambda_z = 0.1
        t_half = nca_half_life(lambda_z)
        expected = math.log(2) / 0.1
        assert abs(t_half - expected) < 1e-10

    def test_auc_0_t(self, init):
        """Test AUC 0-t calculation."""
        auc = auc_0_t(TEST_TIMES, TEST_CONC)
        assert auc > 0.0

    def test_auc_methods(self, init):
        """Test different AUC calculation methods."""
        auc_lin = auc_0_t(TEST_TIMES, TEST_CONC, method="linear")
        auc_log = auc_0_t(TEST_TIMES, TEST_CONC, method="log_linear")
        auc_mixed = auc_0_t(TEST_TIMES, TEST_CONC, method="lin_log_mixed")

        # All should be positive
        assert auc_lin > 0.0
        assert auc_log > 0.0
        assert auc_mixed > 0.0

        # Results should be similar (within 20%)
        assert 0.8 < auc_lin / auc_log < 1.2


# ============================================================================
# PK Parameters Tests
# ============================================================================

class TestPKParameters:
    """Tests for secondary PK parameters."""

    def test_nca_mrt(self, init):
        """Test MRT calculation."""
        mrt = nca_mrt(100.0, 20.0)  # AUMC=100, AUC=20
        assert mrt == 5.0

    def test_nca_mrt_iv_infusion(self, init):
        """Test MRT adjustment for IV infusion."""
        mrt = nca_mrt(100.0, 20.0, route="iv_infusion", t_inf=1.0)
        assert mrt == 4.5  # 5.0 - 1.0/2

    def test_nca_cl_f(self, init):
        """Test CL/F calculation."""
        cl_f = nca_cl_f(100.0, 20.0)  # dose=100, auc=20
        assert cl_f == 5.0

    def test_nca_vz_f(self, init):
        """Test Vz/F calculation."""
        vz_f = nca_vz_f(100.0, 0.1, 20.0)  # dose=100, lambda_z=0.1, auc=20
        assert vz_f == 50.0


# ============================================================================
# Multiple Dose Metrics Tests
# ============================================================================

class TestMultipleDoseMetrics:
    """Tests for multiple dose metrics."""

    def test_nca_accumulation_index(self, init):
        """Test accumulation index calculation."""
        rac = nca_accumulation_index(25.0, 20.0)  # AUC_ss=25, AUC_sd=20
        assert rac == 1.25

    def test_nca_ptf(self, init):
        """Test PTF calculation."""
        ptf = nca_ptf(4.2, 2.5, 3.0)  # Cmax=4.2, Cmin=2.5, Cavg=3.0
        expected = 100.0 * (4.2 - 2.5) / 3.0
        assert abs(ptf - expected) < 1e-10

    def test_nca_swing(self, init):
        """Test Swing calculation."""
        swing = nca_swing(4.2, 2.5)  # Cmax=4.2, Cmin=2.5
        expected = 100.0 * (4.2 - 2.5) / 2.5
        assert abs(swing - expected) < 1e-10


# ============================================================================
# Bioequivalence Tests
# ============================================================================

class TestBioequivalence:
    """Tests for bioequivalence analysis."""

    # Paired crossover data
    test_values = [20.0, 22.0, 18.0, 25.0, 21.0, 19.0, 23.0, 20.0]
    ref_values = [19.0, 21.0, 17.0, 24.0, 20.0, 18.0, 22.0, 19.0]

    def test_geometric_mean(self, init):
        """Test geometric mean calculation."""
        gm = geometric_mean(self.test_values)
        assert gm > 0.0

        import math
        expected = math.exp(sum(math.log(v) for v in self.test_values) / len(self.test_values))
        assert abs(gm - expected) < 1e-10

    def test_geometric_mean_ratio(self, init):
        """Test GMR calculation."""
        gmr = geometric_mean_ratio(self.test_values, self.ref_values)
        assert gmr > 0.0
        assert abs(gmr - 1.05) < 0.1  # Close to 1 for similar values

    def test_within_subject_cv(self, init):
        """Test within-subject CV calculation."""
        cv = within_subject_cv(self.test_values, self.ref_values)
        assert cv > 0.0
        assert cv < 100.0

    def test_bioequivalence_90ci(self, init):
        """Test 90% CI calculation."""
        result = bioequivalence_90ci(self.test_values, self.ref_values)

        assert "gmr" in result
        assert "ci_lower" in result
        assert "ci_upper" in result
        assert "cv_intra" in result
        assert "n" in result

        assert result["gmr"] > 0.0
        assert result["ci_lower"] < result["gmr"]
        assert result["ci_upper"] > result["gmr"]
        assert result["n"] == len(self.test_values)

    def test_tost_analysis(self, init):
        """Test TOST analysis."""
        result = tost_analysis(self.test_values, self.ref_values)

        assert "conclusion" in result
        assert result["conclusion"] in ["bioequivalent", "not_bioequivalent"]

    def test_be_conclusion(self, init):
        """Test BE conclusion determination."""
        # Within limits
        assert be_conclusion(0.85, 1.20) == "bioequivalent"

        # Outside limits
        assert be_conclusion(0.75, 1.10) == "inconclusive"
        assert be_conclusion(0.90, 1.30) == "inconclusive"

        # Completely outside
        assert be_conclusion(0.70, 0.78) == "not_bioequivalent"
        assert be_conclusion(1.30, 1.50) == "not_bioequivalent"


# ============================================================================
# Full NCA Workflow Tests
# ============================================================================

class TestFullNCAWorkflow:
    """Tests for full NCA workflow."""

    def test_run_nca_single_dose(self, init):
        """Test full NCA for single dose."""
        result = run_nca(TEST_TIMES, TEST_CONC, TEST_DOSE)

        assert isinstance(result, NCAResult)
        assert result.cmax == 1.62
        assert result.tmax == 2.0
        assert result.clast == 0.08
        assert result.tlast == 24.0
        assert result.auc_0_t > 0.0

    def test_run_nca_multiple_dose(self, init):
        """Test full NCA for multiple dose."""
        result = run_nca(
            SS_TIMES, SS_CONC, TEST_DOSE,
            dosing_type="multiple",
            tau=SS_TAU
        )

        assert result.cmin is not None
        assert result.cavg is not None
        assert result.auc_0_tau is not None

    def test_run_nca_iv_bolus(self, init):
        """Test full NCA for IV bolus."""
        result = run_nca(
            TEST_TIMES, TEST_CONC, TEST_DOSE,
            route="iv_bolus"
        )

        assert result.cmax > 0.0
        assert result.metadata["route"] == "iv_bolus"

    def test_nca_config(self, init):
        """Test NCA with custom configuration."""
        config = NCAConfig(
            method="log_linear",
            lambda_z_min_points=4,
            lambda_z_r2_threshold=0.85,
        )

        result = run_nca(TEST_TIMES, TEST_CONC, TEST_DOSE, config=config)
        assert result.auc_0_t > 0.0

    def test_nca_result_to_dict(self, init):
        """Test NCAResult to dict conversion."""
        result = run_nca(TEST_TIMES, TEST_CONC, TEST_DOSE)
        result_dict = result.to_dict()

        assert isinstance(result_dict, dict)
        assert "cmax" in result_dict
        assert "auc_0_t" in result_dict
        assert "quality_flags" in result_dict

    def test_quality_flags(self, init):
        """Test that quality flags are properly returned."""
        result = run_nca(TEST_TIMES, TEST_CONC, TEST_DOSE)

        assert isinstance(result.quality_flags, list)
        assert isinstance(result.warnings, list)
        assert isinstance(result.metadata, dict)


# ============================================================================
# Edge Cases Tests
# ============================================================================

class TestEdgeCases:
    """Tests for edge cases."""

    def test_minimum_data_points(self, init):
        """Test with minimum data points."""
        t = [0.0, 1.0, 2.0]
        c = [0.0, 1.0, 0.5]

        result = run_nca(t, c, 100.0)
        assert result.cmax == 1.0

    def test_constant_concentration(self, init):
        """Test with constant concentration."""
        t = [0.0, 1.0, 2.0, 3.0, 4.0]
        c = [1.0, 1.0, 1.0, 1.0, 1.0]

        result = run_nca(t, c, 100.0)
        assert result.cmax == 1.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
