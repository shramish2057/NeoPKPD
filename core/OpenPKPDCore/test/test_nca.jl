# NCA (Non-Compartmental Analysis) Tests
# Tests for FDA/EMA compliant NCA metrics

using Test
using OpenPKPDCore

# =============================================================================
# Test Data Setup
# =============================================================================

# Standard oral PK profile (1-compartment oral first-order absorption)
const TEST_TIMES = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0]
const TEST_CONC = [0.0, 0.82, 1.44, 1.62, 1.28, 0.94, 0.68, 0.36, 0.08]
const TEST_DOSE = 100.0

# IV bolus profile (rapid distribution)
const IV_TIMES = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0]
const IV_CONC = [10.0, 8.2, 6.7, 4.5, 2.0, 0.4, 0.08, 0.016, 0.0006]
const IV_DOSE = 100.0

# Multiple dose steady-state profile
const SS_TIMES = [0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0]
const SS_CONC = [2.5, 4.2, 3.8, 3.0, 2.8, 2.6, 2.5, 2.5]
const SS_TAU = 12.0

# =============================================================================
# NCA Configuration Tests
# =============================================================================

@testset "NCA Configuration" begin
    # Default config - uses LinLogMixedMethod by default
    config = NCAConfig()
    @test config.method isa LinLogMixedMethod
    @test config.lambda_z_min_points == 3
    @test config.lambda_z_r2_threshold == 0.9
    @test config.extrapolation_max_pct == 20.0
    @test config.significant_digits == 3
    @test config.blq_handling isa BLQZero

    # Custom config
    config_custom = NCAConfig(
        method = LogLinearMethod(),
        lambda_z_min_points = 4,
        lambda_z_r2_threshold = 0.95,
        extrapolation_max_pct = 10.0,
        significant_digits = 4,
        blq_handling = BLQLLOQHalf(),
        lloq = 0.1
    )
    @test config_custom.method isa LogLinearMethod
    @test config_custom.lambda_z_min_points == 4
    @test config_custom.lloq == 0.1
end

# =============================================================================
# Exposure Metrics Tests
# =============================================================================

@testset "Exposure Metrics" begin
    @testset "Cmax and Tmax" begin
        cmax = nca_cmax(TEST_CONC)
        tmax = nca_tmax(TEST_TIMES, TEST_CONC)

        @test cmax == 1.62
        @test tmax == 2.0

        # find_cmax returns both value and index
        cmax_val, cmax_idx = find_cmax(TEST_TIMES, TEST_CONC)
        @test cmax_val == 1.62
        @test cmax_idx == 4
    end

    @testset "Cmin" begin
        cmin = nca_cmin(TEST_CONC)
        @test cmin == 0.0  # First point is zero

        # Non-zero minimum
        cmin_ss = nca_cmin(SS_CONC)
        @test cmin_ss == 2.5
    end

    @testset "Clast and Tlast" begin
        clast = nca_clast(TEST_TIMES, TEST_CONC)
        tlast = nca_tlast(TEST_TIMES, TEST_CONC)

        @test clast == 0.08
        @test tlast == 24.0

        # With LLOQ filtering
        clast_lloq = nca_clast(TEST_TIMES, TEST_CONC; lloq=0.1)
        tlast_lloq = nca_tlast(TEST_TIMES, TEST_CONC; lloq=0.1)
        @test clast_lloq == 0.36
        @test tlast_lloq == 12.0
    end

    @testset "find_clast" begin
        clast, tlast, idx = find_clast(TEST_TIMES, TEST_CONC)
        @test clast == 0.08
        @test tlast == 24.0
        @test idx == 9

        # With LLOQ
        clast2, tlast2, idx2 = find_clast(TEST_TIMES, TEST_CONC; lloq=0.5)
        @test clast2 == 0.68
        @test tlast2 == 8.0
        @test idx2 == 7
    end

    @testset "time_above_concentration" begin
        t = [0.0, 1.0, 2.0, 3.0, 4.0]
        c = [0.0, 2.0, 4.0, 2.0, 0.0]

        time_above = time_above_concentration(t, c, 1.0)
        @test time_above ≈ 3.0 atol=0.1  # Approximately 3 hours above threshold
    end
end

# =============================================================================
# Lambda-z Estimation Tests
# =============================================================================

@testset "Lambda-z Estimation" begin
    config = NCAConfig()

    @testset "Basic lambda_z estimation" begin
        result = estimate_lambda_z(TEST_TIMES, TEST_CONC, config)

        @test result.lambda_z !== nothing
        @test result.lambda_z > 0.0
        @test result.t_half !== nothing
        @test result.t_half > 0.0
        @test result.r_squared !== nothing
        @test result.r_squared >= config.lambda_z_r2_threshold
        @test result.n_points >= config.lambda_z_min_points
        @test result.quality_flag in [:good, :warning]
    end

    @testset "Lambda_z relationship to t_half" begin
        result = estimate_lambda_z(TEST_TIMES, TEST_CONC, config)

        if result.lambda_z !== nothing && result.t_half !== nothing
            # t_half = ln(2) / lambda_z
            expected_t_half = log(2) / result.lambda_z
            @test result.t_half ≈ expected_t_half atol=1e-10
        end
    end

    @testset "IV bolus lambda_z" begin
        result = estimate_lambda_z(IV_TIMES, IV_CONC, config)

        @test result.lambda_z !== nothing
        # IV bolus typically has clear terminal phase
        @test result.r_squared !== nothing
        @test result.r_squared > 0.9
    end

    @testset "Insufficient data handling" begin
        # Only 2 points in terminal phase
        short_t = [0.0, 1.0, 2.0]
        short_c = [0.0, 1.0, 0.5]

        result = estimate_lambda_z(short_t, short_c, config)
        # Should handle gracefully
        @test result.quality_flag in [:good, :warning, :insufficient]
    end

    @testset "nca_half_life" begin
        lambda_z = 0.1  # 1/hr
        t_half = nca_half_life(lambda_z)
        @test t_half ≈ log(2) / 0.1 atol=1e-10
    end
end

# =============================================================================
# AUC Calculation Tests
# =============================================================================

@testset "AUC Calculations" begin
    config_linear = NCAConfig(method=LinearMethod())
    config_log = NCAConfig(method=LogLinearMethod())
    config_mixed = NCAConfig(method=LinLogMixedMethod())

    @testset "AUC 0-t (Linear)" begin
        auc = auc_0_t(TEST_TIMES, TEST_CONC, config_linear)
        @test auc > 0.0
        # Manual trapezoid check for first interval
        expected_first = 0.5 * (0.0 + 0.82) * (0.5 - 0.0)
        @test auc > expected_first  # AUC should be larger than first interval
    end

    @testset "AUC 0-t (Log-linear)" begin
        auc = auc_0_t(TEST_TIMES, TEST_CONC, config_log)
        @test auc > 0.0
    end

    @testset "AUC 0-t (Lin-Log Mixed)" begin
        auc = auc_0_t(TEST_TIMES, TEST_CONC, config_mixed)
        @test auc > 0.0
    end

    @testset "AUC methods comparison" begin
        auc_linear = auc_0_t(TEST_TIMES, TEST_CONC, config_linear)
        auc_log = auc_0_t(TEST_TIMES, TEST_CONC, config_log)
        auc_mixed = auc_0_t(TEST_TIMES, TEST_CONC, config_mixed)

        # All methods should give positive results
        @test auc_linear > 0.0
        @test auc_log > 0.0
        @test auc_mixed > 0.0

        # Results should be similar (within 20%)
        @test auc_linear / auc_log > 0.8
        @test auc_linear / auc_log < 1.2
    end

    @testset "AUC 0-inf" begin
        lambda_z_result = estimate_lambda_z(TEST_TIMES, TEST_CONC, config_linear)

        if lambda_z_result.lambda_z !== nothing
            clast, tlast, _ = find_clast(TEST_TIMES, TEST_CONC)
            auc_inf, extra_pct = auc_0_inf(
                TEST_TIMES, TEST_CONC,
                lambda_z_result.lambda_z, clast, config_linear
            )

            @test auc_inf > 0.0
            @test extra_pct >= 0.0
            @test extra_pct < 100.0

            # AUC_inf should be larger than AUC_0-t
            auc_t = auc_0_t(TEST_TIMES, TEST_CONC, config_linear)
            @test auc_inf > auc_t
        end
    end

    @testset "AUC 0-tau" begin
        auc_tau = auc_0_tau(SS_TIMES, SS_CONC, SS_TAU, config_linear)
        @test auc_tau > 0.0

        # Should equal AUC_0-t when tau >= tlast
        auc_t = auc_0_t(SS_TIMES, SS_CONC, config_linear)
        @test auc_tau ≈ auc_t atol=0.01
    end

    @testset "AUMC 0-t" begin
        aumc = aumc_0_t(TEST_TIMES, TEST_CONC, config_linear)
        @test aumc > 0.0
    end

    @testset "AUC partial" begin
        # Partial AUC from 1.0 to 4.0 hours
        partial = auc_partial(TEST_TIMES, TEST_CONC, 1.0, 4.0, config_linear)
        @test partial > 0.0

        # Should be less than total AUC
        total = auc_0_t(TEST_TIMES, TEST_CONC, config_linear)
        @test partial < total
    end
end

# =============================================================================
# PK Parameters Tests
# =============================================================================

@testset "PK Parameters" begin
    config = NCAConfig()

    @testset "MRT calculation" begin
        aumc_inf = 100.0
        auc_inf = 20.0

        mrt = nca_mrt(aumc_inf, auc_inf)
        @test mrt == 5.0

        # IV infusion adjustment
        mrt_inf = nca_mrt(aumc_inf, auc_inf; route=:iv_infusion, t_inf=1.0)
        @test mrt_inf == 4.5  # MRT - t_inf/2
    end

    @testset "CL/F calculation" begin
        cl_f = nca_cl_f(TEST_DOSE, 20.0)  # dose=100, auc=20
        @test cl_f == 5.0
    end

    @testset "CL calculation (IV)" begin
        cl = nca_cl(IV_DOSE, 25.0)
        @test cl == 4.0
    end

    @testset "Vz/F calculation" begin
        vz_f = nca_vz_f(TEST_DOSE, 0.1, 20.0)  # dose=100, lambda_z=0.1, auc=20
        @test vz_f == 50.0
    end

    @testset "Vss calculation" begin
        vss = nca_vss(5.0, 10.0)  # CL=5, MRT=10
        @test vss == 50.0
    end

    @testset "Vss from AUMC" begin
        vss = nca_vss_from_aumc(TEST_DOSE, 20.0, 100.0)  # dose=100, auc=20, aumc=100
        @test vss == 100.0 * 100.0 / 400.0
    end

    @testset "Vc calculation" begin
        vc = nca_vc(IV_DOSE, 10.0)  # dose=100, C0=10
        @test vc == 10.0
    end

    @testset "Bioavailability" begin
        f = nca_bioavailability(20.0, 100.0, 25.0, 100.0)
        @test f ≈ 0.8 atol=1e-10

        # Different doses
        f2 = nca_bioavailability(10.0, 50.0, 25.0, 100.0)
        @test f2 ≈ 0.8 atol=1e-10
    end

    @testset "Mean Absorption Time" begin
        mat = nca_mean_absorption_time(8.0, 5.0)  # MRT_po=8, MRT_iv=5
        @test mat == 3.0
    end
end

# =============================================================================
# Multiple Dose Metrics Tests
# =============================================================================

@testset "Multiple Dose Metrics" begin
    @testset "Accumulation Index" begin
        rac = nca_accumulation_index(25.0, 20.0)  # AUC_ss=25, AUC_sd=20
        @test rac == 1.25
    end

    @testset "Predicted Accumulation" begin
        rac_pred = nca_accumulation_predicted(0.1, 12.0)  # lambda_z=0.1, tau=12
        expected = 1.0 / (1.0 - exp(-0.1 * 12.0))
        @test rac_pred ≈ expected atol=1e-10
    end

    @testset "PTF (Peak-Trough Fluctuation)" begin
        ptf = nca_ptf(4.2, 2.5, 3.0)  # Cmax=4.2, Cmin=2.5, Cavg=3.0
        expected = 100.0 * (4.2 - 2.5) / 3.0
        @test ptf ≈ expected atol=1e-10
    end

    @testset "Swing" begin
        swing = nca_swing(4.2, 2.5)  # Cmax=4.2, Cmin=2.5
        expected = 100.0 * (4.2 - 2.5) / 2.5
        @test swing ≈ expected atol=1e-10
    end

    @testset "Dose Linearity Index" begin
        doses = [25.0, 50.0, 100.0, 200.0]
        aucs = [5.0, 10.0, 20.0, 40.0]  # Perfect linearity

        result = nca_linearity_index(doses, aucs)
        @test result.beta ≈ 1.0 atol=0.01
        @test result.r_squared ≈ 1.0 atol=0.01
        @test result.is_linear == true

        # Non-linear case
        aucs_nonlinear = [5.0, 15.0, 50.0, 180.0]
        result_nl = nca_linearity_index(doses, aucs_nonlinear)
        @test result_nl.beta > 1.0  # Superproportional
    end

    @testset "Time to Steady State" begin
        t_ss = nca_time_to_steady_state(0.1; fraction=0.90)
        expected = -log(1.0 - 0.90) / 0.1
        @test t_ss ≈ expected atol=1e-10

        # 95% steady state
        t_ss_95 = nca_time_to_steady_state(0.1; fraction=0.95)
        @test t_ss_95 > t_ss
    end

    @testset "Doses to Steady State" begin
        n_doses = nca_time_to_steady_state_doses(0.1, 12.0; fraction=0.90)
        @test n_doses >= 1
        @test n_doses isa Int
    end

    @testset "Cavg calculation" begin
        config = NCAConfig()
        cavg = nca_cavg(SS_TIMES, SS_CONC, SS_TAU, config)

        auc_tau = auc_0_tau(SS_TIMES, SS_CONC, SS_TAU, config)
        expected = auc_tau / SS_TAU
        @test cavg ≈ expected atol=1e-10
    end
end

# =============================================================================
# Bioequivalence Tests
# =============================================================================

@testset "Bioequivalence Analysis" begin
    # Test data: paired crossover study
    test_values = [20.0, 22.0, 18.0, 25.0, 21.0, 19.0, 23.0, 20.0]
    ref_values = [19.0, 21.0, 17.0, 24.0, 20.0, 18.0, 22.0, 19.0]

    @testset "Geometric Mean Ratio" begin
        gmr = geometric_mean_ratio(test_values, ref_values)
        @test gmr > 0.0
        @test gmr ≈ 1.05 atol=0.1  # Close to 1 for similar values
    end

    @testset "Geometric Mean" begin
        gm = geometric_mean(test_values)
        expected = exp(sum(log.(test_values)) / length(test_values))
        @test gm ≈ expected atol=1e-10
    end

    @testset "Within-Subject CV" begin
        cv = within_subject_cv(test_values, ref_values)
        @test cv > 0.0
        @test cv < 100.0  # Should be reasonable CV
    end

    @testset "90% Confidence Interval" begin
        result = bioequivalence_90ci(test_values, ref_values)

        @test haskey(result, :gmr)
        @test haskey(result, :ci_lower)
        @test haskey(result, :ci_upper)
        @test haskey(result, :cv_intra)
        @test haskey(result, :n)

        @test result.gmr > 0.0
        @test result.ci_lower < result.gmr
        @test result.ci_upper > result.gmr
        @test result.n == length(test_values)
    end

    @testset "TOST Analysis" begin
        result = tost_analysis(test_values, ref_values)

        @test result.parameter == :generic
        @test result.be_conclusion in [:bioequivalent, :not_bioequivalent]

        # Check test statistics exist
        @test !isnan(result.t_lower)
        @test !isnan(result.t_upper)
    end

    @testset "BE Conclusion" begin
        # Within limits
        @test be_conclusion(0.85, 1.20) == :bioequivalent

        # Outside limits
        @test be_conclusion(0.75, 1.10) == :inconclusive
        @test be_conclusion(0.90, 1.30) == :inconclusive

        # Completely outside
        @test be_conclusion(0.70, 0.78) == :not_bioequivalent
        @test be_conclusion(1.30, 1.50) == :not_bioequivalent
    end

    @testset "Custom BE Limits" begin
        # Highly variable drug (wider limits)
        conclusion = be_conclusion(0.72, 1.35; theta_lower=0.70, theta_upper=1.43)
        @test conclusion == :bioequivalent
    end

    @testset "Create BE Result" begin
        result = create_be_result(:cmax, test_values, ref_values)

        @test result.parameter == :cmax
        @test result.n_test == length(test_values)
        @test result.n_reference == length(ref_values)
        @test result.gmr > 0.0
        @test result.be_conclusion in [:bioequivalent, :not_bioequivalent, :inconclusive]
    end
end

# =============================================================================
# Full NCA Workflow Tests
# =============================================================================

@testset "Full NCA Workflow" begin
    config = NCAConfig()

    @testset "run_nca - Single Dose Oral" begin
        result = run_nca(TEST_TIMES, TEST_CONC, TEST_DOSE; config=config)

        # Primary exposure
        @test result.cmax == 1.62
        @test result.tmax == 2.0
        @test result.clast == 0.08
        @test result.tlast == 24.0

        # AUC
        @test result.auc_0_t > 0.0

        # Terminal phase (if available)
        @test result.lambda_z_result isa LambdaZResult

        # Quality
        @test result.quality_flags isa Vector{Symbol}
        @test result.warnings isa Vector{String}

        # Metadata
        @test result.metadata["dose"] == TEST_DOSE
        @test result.metadata["route"] == "extravascular"
        @test result.metadata["dosing_type"] == "single"
    end

    @testset "run_nca - IV Bolus" begin
        result = run_nca(IV_TIMES, IV_CONC, IV_DOSE;
                        config=config, route=:iv_bolus)

        @test result.cmax > 0.0
        @test result.auc_0_t > 0.0
        @test result.metadata["route"] == "iv_bolus"
    end

    @testset "run_nca - Multiple Dose" begin
        result = run_nca(SS_TIMES, SS_CONC, TEST_DOSE;
                        config=config, dosing_type=:multiple, tau=SS_TAU)

        @test result.cmin !== nothing
        @test result.cavg !== nothing
        @test result.auc_0_tau !== nothing
        @test result.metadata["tau"] == SS_TAU
    end

    @testset "run_nca - Steady State" begin
        result = run_nca(SS_TIMES, SS_CONC, TEST_DOSE;
                        config=config, dosing_type=:steady_state, tau=SS_TAU)

        @test result.cmin !== nothing
        @test result.ptf !== nothing || result.cmin == 0.0
        @test result.swing !== nothing || result.cmin == 0.0
    end

    @testset "run_nca with custom config" begin
        custom_config = NCAConfig(
            method = LogLinearMethod(),
            lambda_z_r2_threshold = 0.85,
            extrapolation_max_pct = 30.0
        )

        result = run_nca(TEST_TIMES, TEST_CONC, TEST_DOSE; config=custom_config)
        @test result.auc_0_t > 0.0
    end

    @testset "nca_from_simresult" begin
        # Create mock SimResult-like structure
        # This would normally come from a simulation
        t = TEST_TIMES
        c = TEST_CONC

        # Just test the individual NCA works with same data
        result = run_nca(t, c, TEST_DOSE)
        @test result.cmax > 0.0
    end
end

# =============================================================================
# Edge Cases and Error Handling
# =============================================================================

@testset "Edge Cases" begin
    config = NCAConfig()

    @testset "Minimum data points" begin
        t = [0.0, 1.0, 2.0]
        c = [0.0, 1.0, 0.5]

        # Should not throw with 3 points
        result = run_nca(t, c, 100.0; config=config)
        @test result.cmax == 1.0
    end

    @testset "Constant concentration" begin
        t = [0.0, 1.0, 2.0, 3.0, 4.0]
        c = [1.0, 1.0, 1.0, 1.0, 1.0]

        result = run_nca(t, c, 100.0; config=config)
        @test result.cmax == 1.0
        # cmin is only calculated for multiple dose, not single dose
        @test result.cmin === nothing
    end

    @testset "All zeros after peak" begin
        t = [0.0, 1.0, 2.0, 3.0, 4.0]
        c = [0.0, 1.0, 0.0, 0.0, 0.0]

        # Lambda_z estimation may fail, but NCA should complete
        result = run_nca(t, c, 100.0; config=config)
        @test result.cmax == 1.0
    end

    @testset "BLQ handling" begin
        config_blq = NCAConfig(lloq=0.5, blq_handling=BLQZero())

        t = [0.0, 1.0, 2.0, 3.0, 4.0]
        c = [0.0, 2.0, 1.0, 0.3, 0.1]  # Last 2 below LLOQ

        result = run_nca(t, c, 100.0; config=config_blq)
        @test result.cmax == 2.0
    end
end

# =============================================================================
# Input Validation Tests
# =============================================================================

@testset "Input Validation" begin
    config = NCAConfig()

    @testset "Mismatched lengths" begin
        @test_throws AssertionError run_nca([0.0, 1.0], [0.0], 100.0)
    end

    @testset "Negative dose" begin
        @test_throws AssertionError run_nca(TEST_TIMES, TEST_CONC, -100.0)
    end

    @testset "Unsorted times" begin
        @test_throws AssertionError run_nca([1.0, 0.0, 2.0], [1.0, 0.0, 0.5], 100.0)
    end

    @testset "Insufficient points" begin
        @test_throws AssertionError run_nca([0.0, 1.0], [0.0, 1.0], 100.0)
    end

    @testset "Multiple dose without tau" begin
        @test_throws AssertionError run_nca(
            TEST_TIMES, TEST_CONC, TEST_DOSE;
            dosing_type=:multiple
        )
    end
end

# =============================================================================
# Numerical Precision Tests
# =============================================================================

@testset "Numerical Precision" begin
    @testset "AUC precision" begin
        # Simple trapezoid - exact answer is 2.0
        t = [0.0, 1.0, 2.0]
        c = [0.0, 2.0, 2.0]
        config = NCAConfig(method=LinearMethod())

        auc = auc_0_t(t, c, config)
        @test auc ≈ 3.0 atol=1e-10  # Triangle (1.0) + Rectangle (2.0)
    end

    @testset "Log-linear extrapolation" begin
        # Verify extrapolation is reasonable
        t = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
        c = [0.0, 1.0, 0.8, 0.64, 0.512, 0.4096, 0.328, 0.262, 0.21]
        config = NCAConfig()

        result = estimate_lambda_z(t, c, config)

        if result.lambda_z !== nothing
            # Should be approximately -log(0.8) ≈ 0.223
            @test result.lambda_z > 0.15
            @test result.lambda_z < 0.30
        end
    end

    @testset "Round NCA result" begin
        @test round_nca_result(1.2345, 3) ≈ 1.23 atol=0.01
        @test round_nca_result(0.001234, 3) ≈ 0.00123 atol=0.00001
        @test round_nca_result(1234.5, 3) ≈ 1230.0 atol=1.0
    end
end
