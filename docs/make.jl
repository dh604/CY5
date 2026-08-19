using Documenter, GW_CY5
using DocumenterCitations

DocMeta.setdocmeta!(GW_CY5, :DocTestSetup, :(using Oscar, GKMtools, GW_CY5); recursive=true)

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:alpha
)

makedocs(
    sitename = "Membranes, Maps and Sheaves",
    modules  = [GW_CY5],
    plugins  = [bib],
    warnonly = true,
    format   = Documenter.HTML(mathengine = MathJax3()),
    pages    = [
        "Home"       => "index.md",
        "Gromov-Witten computations" => [
            "The package GW_CY5.jl" => "cy5_overview.md",
            "Pipeline"   => "pipeline.md",
            "Spaces" => "spaces.md",
            "Examples"   => [
                    "Julia code in the paper" => "code_from_paper.md",
                    "Details of the Pipeline" => "example.md"
                ],
            "Evidence" => [
                "Overview" => "conj/overview.md",
                "Conjecture 4.2" => "conj/Conj_4_2.md",
                "Section 4.3" => "conj/Conj_4_3.md",
                "Conjecture 5.1" => "conj/Conj_5_1.md",
                "Conjecture 6.1" => "conj/Conj_6_1.md",
                "Conjecture 7.1" => "conj/Conj_7_1.md",
                "Section 7.2" => "conj/Conj_7_2.md",
                "Section 7.3" => "conj/Conj_7_3.md",
                "Section 7.4" => "conj/Conj_7_4.md",
                "Section 7.5" => "conj/Conj_7_5.md",
                "Conjecture 8.1" => "conj/Conj_8_1.md"
            ]
        ],
        "Pandharipande-Thomas computations" => [
            "Overview" => "pt/overview.md"
            "Other page" => "pt/other_page.md"
        ],
        "References" => "references.md"
    ];
    doctest=true
)
