using Janus
using Documenter

DocMeta.setdocmeta!(Janus, :DocTestSetup, :(using Janus); recursive=true)

makedocs(;
    modules=[Janus],
    authors="Jarvist Moore Frost <jarvist@gmail.com> and contributors",
    sitename="Janus.jl",
    format=Documenter.HTML(;
        canonical="https://jarvist.github.io/Janus.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jarvist/Janus.jl",
    devbranch="main",
)
