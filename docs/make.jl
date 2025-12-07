using DBF
using Documenter

DocMeta.setdocmeta!(DBF, :DocTestSetup, :(using DBF); recursive=true)

makedocs(;
    modules=[DBF],
    authors="Nick Mayhall and contributors",
    sitename="DBF.jl",
    format=Documenter.HTML(;
        canonical="https://nmayhall.github.io/DBF.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Theory and Mathematics" => "theory.md",
        "DBF Variants Comparison" => "dbf_variants.md",
        "Repository Structure" => "structure.md",
        "User Guide and Examples" => "guide.md",
    ],
)

deploydocs(;
    repo="github.com/nmayhall/DBF.jl",
    devbranch="main",
)
