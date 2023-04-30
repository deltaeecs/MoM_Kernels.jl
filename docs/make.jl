using MoM_Kernels
using Documenter

DocMeta.setdocmeta!(MoM_Kernels, :DocTestSetup, :(using MoM_Kernels); recursive=true)

makedocs(;
    modules=[MoM_Kernels],
    authors="deltaeecs <1225385871@qq.com> and contributors",
    repo="https://github.com/deltaeecs/MoM_Kernels.jl/blob/{commit}{path}#{line}",
    sitename="MoM_Kernels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://deltaeecs.github.io/MoM_Kernels.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/deltaeecs/MoM_Kernels.jl",
)
