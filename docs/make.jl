using Documenter, CY5
using DocumenterCitations

DocMeta.setdocmeta!(CY5, :DocTestSetup, :(using Oscar, GKMtools, CY5); recursive=true)

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:alpha
)

makedocs(
    sitename = "CY5",
    modules  = [CY5],
    plugins  = [bib],
    warnonly = true,
    format   = Documenter.HTML(mathengine = MathJax3()),
    pages    = [
        "Home"       => "index.md",
        "Gromov-Witten computations" => [
            "Pipeline"   => "pipeline.md",
            "Spaces" => "spaces.md",
            "Example"   => "example.md",
            "Conjectures" => [
                "Overview of Conjectures" => "conj/overview.md",
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
    doctest=false
)
