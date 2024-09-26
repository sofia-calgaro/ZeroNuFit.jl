using Pkg
Pkg.activate("docs")
Pkg.instantiate()

using Documenter
using ZeroNuFit  # Replace with your package's name

makedocs(
    sitename = "ZeroNuFit.jl",
    authors = "S. Calgaro, T. Dixon",
    format = :html,
    pages = [
        "Home" => "index.md",
    ]
)