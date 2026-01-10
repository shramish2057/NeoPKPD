export ARTIFACT_SCHEMA_VERSION

# Schema version 1.1.0 adds compliance_metadata block for FDA 21 CFR Part 11 alignment
# - audit_record: execution ID, timestamps, user/host info
# - integrity: SHA-256 hashes for content, inputs, outputs
# - environment: complete environment snapshot
const ARTIFACT_SCHEMA_VERSION = "1.1.0"

# Schema version history:
# 1.0.0 - Initial release
# 1.1.0 - Added compliance_metadata block for regulatory compliance
