using NeoPKPDCore

root = normpath(joinpath(@__DIR__, "..", ".."))
version_file = joinpath(root, "VERSION")

v = strip(read(version_file, String))

if v != NeoPKPDCore.NEOPKPD_VERSION
    error("VERSION file ($(v)) does not match NeoPKPDCore.NEOPKPD_VERSION ($(NeoPKPDCore.NEOPKPD_VERSION))")
end

println("Version consistency check passed: " * v)
