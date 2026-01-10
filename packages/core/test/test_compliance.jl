using Test
using NeoPKPDCore
using JSON

@testset "FDA 21 CFR Part 11 Compliance Module" begin

    @testset "ComplianceConfig" begin
        @testset "Default configuration" begin
            config = get_compliance_config()
            @test config.level == COMPLIANCE_STANDARD
            @test config.enable_audit_trail == true
            @test config.enable_hashing == true
            @test config.enable_environment_capture == true
            @test config.enable_validation_reports == false
        end

        @testset "Configuration from level" begin
            disabled = ComplianceConfig(COMPLIANCE_DISABLED)
            @test disabled.enable_audit_trail == false
            @test disabled.enable_hashing == false

            minimal = ComplianceConfig(COMPLIANCE_MINIMAL)
            @test minimal.enable_environment_capture == true
            @test minimal.enable_hashing == false

            standard = ComplianceConfig(COMPLIANCE_STANDARD)
            @test standard.enable_audit_trail == true
            @test standard.enable_hashing == true

            strict = ComplianceConfig(COMPLIANCE_STRICT)
            @test strict.enable_validation_reports == true
        end

        @testset "Set and get configuration" begin
            old_config = get_compliance_config()

            set_compliance_config!(COMPLIANCE_DISABLED)
            @test get_compliance_config().level == COMPLIANCE_DISABLED

            set_compliance_config!(COMPLIANCE_STANDARD)
            @test get_compliance_config().level == COMPLIANCE_STANDARD

            # Restore original
            set_compliance_config!(old_config)
        end

        @testset "with_compliance context" begin
            @test get_compliance_config().level == COMPLIANCE_STANDARD

            with_compliance(COMPLIANCE_DISABLED) do
                @test get_compliance_config().level == COMPLIANCE_DISABLED
            end

            # Should be restored
            @test get_compliance_config().level == COMPLIANCE_STANDARD
        end

        @testset "Helper functions" begin
            set_compliance_config!(COMPLIANCE_STANDARD)
            @test is_audit_trail_enabled() == true
            @test is_hashing_enabled() == true
            @test is_environment_capture_enabled() == true
            @test is_validation_reports_enabled() == false

            with_compliance(COMPLIANCE_DISABLED) do
                @test is_audit_trail_enabled() == false
                @test is_hashing_enabled() == false
            end
        end
    end

    @testset "EnvironmentSnapshot" begin
        @testset "Capture environment" begin
            env = capture_environment()

            @test !isempty(env.julia_version)
            @test !isempty(env.neopkpd_version)
            @test env.neopkpd_version == NEOPKPD_VERSION
            @test env.artifact_schema_version == ARTIFACT_SCHEMA_VERSION
            @test !isempty(env.os)
            @test env.cpu_cores > 0
            @test env.memory_gb > 0
            @test !isempty(env.architecture)
            @test env.word_size in [32, 64]
            @test !isempty(env.capture_timestamp)
        end

        @testset "Serialize environment" begin
            env = capture_environment()
            serialized = serialize_environment(env)

            @test serialized["julia_version"] == env.julia_version
            @test serialized["neopkpd_version"] == env.neopkpd_version
            @test serialized["cpu_cores"] == env.cpu_cores
        end

        @testset "Deserialize environment" begin
            env = capture_environment()
            serialized = serialize_environment(env)
            deserialized = deserialize_environment(serialized)

            @test deserialized.julia_version == env.julia_version
            @test deserialized.neopkpd_version == env.neopkpd_version
            @test deserialized.cpu_cores == env.cpu_cores
        end
    end

    @testset "IntegrityMetadata" begin
        @testset "Compute content hash" begin
            hash1 = compute_content_hash("test data")
            hash2 = compute_content_hash("test data")
            hash3 = compute_content_hash("different data")

            @test length(hash1) == 64
            @test all(c -> c in "0123456789abcdef", hash1)
            @test hash1 == hash2
            @test hash1 != hash3
        end

        @testset "Compute dict hash" begin
            d = Dict("key" => "value", "number" => 42)
            hash = compute_content_hash(d)

            @test length(hash) == 64
        end

        @testset "Canonical JSON determinism" begin
            d1 = Dict("b" => 2, "a" => 1)
            d2 = Dict("a" => 1, "b" => 2)

            hash1 = compute_content_hash(d1)
            hash2 = compute_content_hash(d2)

            @test hash1 == hash2  # Should be same despite insertion order
        end

        @testset "Compute integrity metadata" begin
            artifact = Dict{String,Any}(
                "model_spec" => Dict("CL" => 1.0, "V" => 10.0),
                "result" => Dict("concentration" => [10.0, 5.0, 2.5]),
                "semantics_fingerprint" => Dict("version" => "1.0.0")
            )

            meta = compute_integrity_metadata(artifact)

            @test length(meta.content_hash) == 64
            @test meta.hash_algorithm == "SHA-256"
            @test length(meta.input_hash) == 64
            @test length(meta.output_hash) == 64
            @test length(meta.semantics_hash) == 64
            @test !isempty(meta.hash_computed_at)
        end

        @testset "Serialize integrity" begin
            artifact = Dict{String,Any}("key" => "value")
            meta = compute_integrity_metadata(artifact)
            serialized = serialize_integrity(meta)

            @test serialized["content_hash"] == meta.content_hash
            @test serialized["hash_algorithm"] == "SHA-256"
        end
    end

    @testset "AuditRecord" begin
        @testset "Create audit record" begin
            record = create_audit_record()

            @test !isempty(record.execution_id)
            @test length(split(record.execution_id, "-")) == 5  # UUID format
            @test !isempty(record.timestamp_utc)
            @test endswith(record.timestamp_utc, "Z")  # UTC
            @test !isempty(record.system_user)
            @test !isempty(record.hostname)
            @test record.action == AUDIT_CREATE
            @test record.previous_execution_id === nothing
            @test record.neopkpd_version == NEOPKPD_VERSION
            @test length(record.record_hash) == 64
        end

        @testset "Create replay audit record" begin
            original = create_audit_record()
            replay = create_audit_record(
                action=AUDIT_REPLAY,
                previous_execution_id=original.execution_id
            )

            @test replay.action == AUDIT_REPLAY
            @test replay.previous_execution_id == original.execution_id
            @test replay.execution_id != original.execution_id
        end

        @testset "Verify audit record" begin
            record = create_audit_record()
            @test verify_audit_record(record) == true

            # Tampered record should fail
            tampered = AuditRecord(
                record.execution_id,
                record.timestamp_utc,
                record.timezone_offset,
                "fake_user",  # Tampered
                record.hostname,
                record.execution_context,
                record.action,
                record.previous_execution_id,
                record.neopkpd_version,
                record.record_hash  # Original hash won't match
            )
            @test verify_audit_record(tampered) == false
        end

        @testset "Serialize audit record" begin
            record = create_audit_record()
            serialized = serialize_audit_record(record)

            @test serialized["execution_id"] == record.execution_id
            @test serialized["action"] == "create"
            @test serialized["neopkpd_version"] == NEOPKPD_VERSION
        end

        @testset "Deserialize audit record" begin
            record = create_audit_record()
            serialized = serialize_audit_record(record)
            deserialized = deserialize_audit_record(serialized)

            @test deserialized.execution_id == record.execution_id
            @test deserialized.action == record.action
        end
    end

    @testset "Compliance Metadata Integration" begin
        @testset "Add compliance metadata to artifact" begin
            spec = ModelSpec(
                OneCompIVBolus(),
                "test",
                OneCompIVBolusParams(1.0, 10.0),
                [DoseEvent(0.0, 100.0)]
            )
            grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
            solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

            result = simulate(spec, grid, solver)
            artifact = serialize_execution(
                model_spec=spec,
                grid=grid,
                solver=solver,
                result=result
            )

            add_compliance_metadata!(artifact)

            @test haskey(artifact, "compliance_metadata")
            meta = artifact["compliance_metadata"]
            @test haskey(meta, "audit_record")
            @test haskey(meta, "integrity")
            @test haskey(meta, "environment")
        end

        @testset "Verify artifact integrity" begin
            spec = ModelSpec(
                OneCompIVBolus(),
                "test",
                OneCompIVBolusParams(1.0, 10.0),
                [DoseEvent(0.0, 100.0)]
            )
            grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
            solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

            result = simulate(spec, grid, solver)
            artifact = serialize_execution(
                model_spec=spec,
                grid=grid,
                solver=solver,
                result=result
            )
            add_compliance_metadata!(artifact)

            # Should be valid
            verification = verify_artifact_integrity(artifact)
            @test verification.is_valid == true
            @test verification.content_hash_valid == true

            # Tamper and verify (should fail)
            artifact["result"]["observations"]["conc"][1] = 999.0
            verification2 = verify_artifact_integrity(artifact)
            @test verification2.is_valid == false
        end

        @testset "Has compliance metadata" begin
            artifact_with = Dict{String,Any}("compliance_metadata" => Dict())
            artifact_without = Dict{String,Any}("data" => 123)

            @test has_compliance_metadata(artifact_with) == true
            @test has_compliance_metadata(artifact_without) == false
        end

        @testset "Get execution ID" begin
            artifact = Dict{String,Any}(
                "compliance_metadata" => Dict{String,Any}(
                    "audit_record" => Dict{String,Any}(
                        "execution_id" => "test-id-123"
                    )
                )
            )

            @test get_execution_id(artifact) == "test-id-123"
            @test get_execution_id(Dict{String,Any}()) === nothing
        end

        @testset "Compliance disabled mode" begin
            artifact = Dict{String,Any}("data" => 123)

            with_compliance(COMPLIANCE_DISABLED) do
                add_compliance_metadata!(artifact)
            end

            @test !haskey(artifact, "compliance_metadata")
        end
    end

    @testset "Validation Reports" begin
        @testset "Generate validation report" begin
            tests = [
                ValidationTestResult(
                    "TEST-001", "Test 1", "Description 1",
                    VALIDATION_PASSED, "expected", "actual", nothing, 100.0
                ),
                ValidationTestResult(
                    "TEST-002", "Test 2", "Description 2",
                    VALIDATION_FAILED, "expected", "wrong", "Error message", 50.0
                ),
            ]

            report = generate_validation_report(
                report_type="IQ",
                tests=tests,
                notes=["Test note"],
                execution_time_ms=150.0
            )

            @test report.report_type == "IQ"
            @test report.total_tests == 2
            @test report.passed_tests == 1
            @test report.failed_tests == 1
            @test report.overall_status == VALIDATION_FAILED
        end

        @testset "Write validation report" begin
            tests = [
                ValidationTestResult(
                    "TEST-001", "Test 1", "Description",
                    VALIDATION_PASSED, "expected", "actual", nothing, 100.0
                ),
            ]

            report = generate_validation_report(
                report_type="IQ",
                tests=tests
            )

            temp_path = tempname() * ".md"
            write_validation_report_markdown(report, temp_path)

            @test isfile(temp_path)
            content = read(temp_path, String)
            @test occursin("IQ Validation Report", content)
            @test occursin("TEST-001", content)

            rm(temp_path)
        end
    end

    @testset "Schema Version" begin
        @test ARTIFACT_SCHEMA_VERSION == "1.1.0"
    end

end
