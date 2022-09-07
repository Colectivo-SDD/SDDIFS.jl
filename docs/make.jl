push!(LOAD_PATH, "./../src/")
using SDDIFS

using Documenter
makedocs(
    modules = [SDDIFS],
    sitename = "SDDIFS Reference",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        warn_outdated = true,
        collapselevel=1,
        )
)
