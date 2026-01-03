@testset "Trial Module" begin

    @testset "Study Designs" begin
        @testset "Parallel Design" begin
            # Basic parallel design
            design = parallel_design(2)
            @test design isa ParallelDesign
            @test design.n_arms == 2
            @test length(design.randomization_ratio) == 2
            @test sum(design.randomization_ratio) ≈ 1.0

            # 3-arm with custom ratio
            design3 = parallel_design(3, randomization_ratio=[0.5, 0.25, 0.25])
            @test design3.n_arms == 3
            @test design3.randomization_ratio ≈ [0.5, 0.25, 0.25]
        end

        @testset "Crossover Design" begin
            # 2x2 crossover
            design = crossover_2x2()
            @test design isa CrossoverDesign
            @test design.n_periods == 2
            @test design.n_sequences == 2
            @test design.washout_duration == 7.0

            # Williams design
            williams = williams_design(4)
            @test williams.n_periods == 4
            @test length(williams.sequence_assignments) == 8
        end

        @testset "Dose Escalation" begin
            # 3+3 design
            dose_levels = [10.0, 25.0, 50.0, 100.0, 200.0]
            design = dose_escalation_3plus3(dose_levels)
            @test design isa DoseEscalationDesign
            @test design.starting_dose == 10.0
            @test design.cohort_size == 3
            @test design.escalation_rule isa ThreePlusThree

            # mTPI design
            mtpi = dose_escalation_mtpi(dose_levels, target_dlt_rate=0.30)
            @test mtpi.escalation_rule isa mTPI
            @test mtpi.escalation_rule.target_dlt_rate == 0.30

            # CRM design
            crm = dose_escalation_crm(dose_levels)
            @test crm.escalation_rule isa CRM
        end

        @testset "Bioequivalence Design" begin
            design = bioequivalence_design()
            @test design isa BioequivalenceDesign
            @test design.bioequivalence_limits == (0.80, 1.25)
            @test :cmax in design.parameters
            @test :auc_0_inf in design.parameters
        end

        @testset "Adaptive Design" begin
            base = parallel_design(2)
            adaptive = adaptive_design(base, interim_analyses=[0.5])
            @test adaptive isa AdaptiveDesign
            @test adaptive.interim_analyses == [0.5]
            @test adaptive.alpha_spending == :obrien_fleming
        end

        @testset "Design Descriptions" begin
            @test contains(get_design_description(parallel_design(2)), "parallel")
            @test contains(get_design_description(crossover_2x2()), "crossover")
            @test contains(get_design_description(bioequivalence_design()), "Bioequivalence")
        end
    end

    @testset "Dosing Regimens" begin
        @testset "Standard Regimens" begin
            qd = dosing_qd(100.0, 7)
            @test qd.frequency isa QD
            @test qd.dose_amount == 100.0
            @test qd.duration_days == 7

            bid = dosing_bid(50.0, 14)
            @test bid.frequency isa BID
            @test length(bid.dose_times) == 2

            tid = dosing_tid(25.0, 7)
            @test tid.frequency isa TID
            @test length(tid.dose_times) == 3

            # With loading dose
            with_loading = dosing_qd(100.0, 7, loading_dose=200.0)
            @test with_loading.loading_dose == 200.0
        end

        @testset "Dose Event Times" begin
            regimen = dosing_bid(50.0, 3)
            times = dose_event_times(regimen)
            @test length(times) == 6  # 3 days * 2 doses/day
            @test times[1] == 8.0  # First dose at 8 AM
        end

        @testset "Titration Regimen" begin
            tit = titration_regimen(25.0, 100.0, 4, 7)
            @test tit isa TitrationRegimen
            @test length(tit.steps) == 4
            @test tit.steps[1].dose == 25.0
            @test tit.steps[end].dose == 100.0

            duration = total_regimen_duration(tit)
            @test duration == 28  # 4 steps * 7 days
        end

        @testset "Generate Doses" begin
            regimen = dosing_qd(100.0, 5)
            doses = generate_doses(regimen)
            @test length(doses) == 5
            @test all(d == 100.0 for d in doses)

            # With loading dose
            regimen_load = dosing_qd(100.0, 5, loading_dose=200.0)
            doses_load = generate_doses(regimen_load)
            @test doses_load[1] == 200.0
            @test doses_load[2] == 100.0
        end
    end

    @testset "Virtual Population" begin
        @testset "Demographic Specs" begin
            demo = default_demographic_spec()
            @test demo.age_mean == 35.0
            @test demo.weight_mean == 75.0

            hv = healthy_volunteer_spec()
            @test hv.age_range == (18.0, 45.0)
        end

        @testset "Patient Populations" begin
            demo, disease = patient_population_spec(:diabetes)
            @test demo.age_mean > 40.0  # Older population
            @test disease.name == :diabetes

            demo2, disease2 = patient_population_spec(:renal)
            @test disease2.name == :renal_impairment
        end

        @testset "Population Generation" begin
            spec = VirtualPopulationSpec(
                demographics = DemographicSpec(age_mean=40.0, age_sd=10.0),
                seed = UInt64(42)
            )
            pop = generate_virtual_population(spec, 50)

            @test length(pop) == 50
            @test all(p.id >= 1 for p in pop)
            @test all(18.0 <= p.age <= 75.0 for p in pop)
            @test all(p.sex in [:male, :female] for p in pop)
        end

        @testset "Population Summary" begin
            spec = VirtualPopulationSpec(seed = UInt64(42))
            pop = generate_virtual_population(spec, 100)
            summary = summarize_population(pop)

            @test summary[:n] == 100
            @test haskey(summary, :age_mean)
            @test haskey(summary, :weight_mean)
            @test haskey(summary, :female_proportion)
            @test 0.0 <= summary[:female_proportion] <= 1.0
        end
    end

    @testset "Trial Events" begin
        @testset "Dropout Simulation" begin
            spec = DropoutSpec(random_rate_per_day=0.01)
            dropouts = simulate_dropout(spec, 30.0, 100)

            @test length(dropouts) <= 100
            @test all(d.time_days <= 30.0 for d in dropouts)
            @test all(d.reason in [:random, :ae, :non_compliance] for d in dropouts)
        end

        @testset "Compliance" begin
            spec = ComplianceSpec(mean_compliance=0.85)
            compliance = apply_compliance(spec, 50, 28.0)

            @test length(compliance) == 50
            @test all(0.0 <= c <= 1.0 for c in compliance)
            @test 0.7 < sum(compliance)/50 < 1.0  # Average around 0.85
        end

        @testset "Survival Time" begin
            dropouts = [DropoutEvent(1, 10.0, :random), DropoutEvent(3, 15.0, :ae)]

            time1, comp1 = calculate_survival_time(dropouts, 1, 30.0)
            @test time1 == 10.0
            @test !comp1

            time2, comp2 = calculate_survival_time(dropouts, 2, 30.0)
            @test time2 == 30.0
            @test comp2
        end
    end

    @testset "Endpoints" begin
        @testset "PK Endpoint Extraction" begin
            result = Dict{String, Any}(
                "t" => [0.0, 1.0, 2.0, 4.0, 8.0],
                "observations" => Dict{String, Vector{Float64}}(
                    "conc" => [0.0, 10.0, 8.0, 5.0, 2.0]
                )
            )

            endpoint_cmax = PKEndpoint(:cmax, metric=:cmax)
            cmax = extract_pk_endpoint(result, endpoint_cmax)
            @test cmax == 10.0

            endpoint_tmax = PKEndpoint(:tmax, metric=:tmax)
            tmax = extract_pk_endpoint(result, endpoint_tmax)
            @test tmax == 1.0

            endpoint_auc = PKEndpoint(:auc, metric=:auc_0_t)
            auc = extract_pk_endpoint(result, endpoint_auc)
            @test auc > 0
        end

        @testset "Endpoint Analysis" begin
            values = [10.0, 12.0, 9.0, 11.0, 13.0, 8.0, 10.5, 11.5, 9.5, 10.0]
            endpoint = PKEndpoint(:test, metric=:cmax)
            stats = analyze_endpoint(values, endpoint)

            @test stats[:n] == 10
            @test 9.0 < stats[:mean] < 12.0
            @test stats[:min] == 8.0
            @test stats[:max] == 13.0
        end

        @testset "Arm Comparison" begin
            arm1 = [10.0, 11.0, 12.0, 9.0, 10.5, 11.5, 10.0, 9.5]
            arm2 = [14.0, 15.0, 16.0, 13.0, 14.5, 15.5, 14.0, 13.5]

            comparison = compare_arms(arm1, arm2, test=:ttest)
            @test comparison[:difference] > 0  # arm2 > arm1
            @test haskey(comparison, :t_statistic)
            @test haskey(comparison, :ci_lower)
            @test haskey(comparison, :ci_upper)
        end

        @testset "Responder Analysis" begin
            values = [0.3, 0.6, 0.8, 0.4, 0.9, 0.2, 0.7, 0.5, 0.85, 0.35]
            result = responder_analysis(values, 0.5, direction=:greater)

            @test result[:n] == 10
            @test result[:n_responders] == 6  # values >= 0.5
            @test result[:response_rate] == 0.6
        end
    end

    @testset "Power Analysis" begin
        @testset "Analytical Power" begin
            # Standard two-sample t-test scenario
            power = estimate_power_analytical(50, 0.5, 1.0)
            @test 0.0 < power < 1.0

            # Larger effect = higher power
            power_large = estimate_power_analytical(50, 1.0, 1.0)
            @test power_large > power

            # Larger N = higher power
            power_large_n = estimate_power_analytical(100, 0.5, 1.0)
            @test power_large_n > power
        end

        @testset "Sample Size Estimation" begin
            result = estimate_sample_size(0.80, 0.5, 1.0)
            @test result isa SampleSizeResult
            @test result.n_per_arm > 0
            @test result.achieved_power >= result.target_power
        end

        @testset "Alpha Spending" begin
            # O'Brien-Fleming - conservative early
            alpha_50 = alpha_spending_function(0.5, 0.05, :obrien_fleming)
            @test alpha_50 < 0.025  # Less than half of total alpha

            # Linear spending
            alpha_lin = alpha_spending_function(0.5, 0.05, :linear)
            @test alpha_lin ≈ 0.025 atol=0.001

            # At final analysis
            alpha_final = alpha_spending_function(1.0, 0.05, :obrien_fleming)
            @test alpha_final ≈ 0.05 atol=0.01
        end

        @testset "Incremental Alpha" begin
            fractions = [0.5, 1.0]
            incremental = incremental_alpha(fractions, 0.05, :obrien_fleming)

            @test length(incremental) == 2
            @test sum(incremental) ≈ 0.05 atol=0.01
        end
    end

end
