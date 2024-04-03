using OptQit
using Documenter

DocMeta.setdocmeta!(OptQit, :DocTestSetup, :(using OptQit); recursive=true)

makedocs(;
    modules=[OptQit],
    authors="Yusheng Zhao <yushengzhao2020@outlook.com> and contributors",
    sitename="OptQit.jl",
    format=Documenter.HTML(;
        canonical="https://exAClior.github.io/OptQit.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/exAClior/OptQit.jl",
    devbranch="main",
)
