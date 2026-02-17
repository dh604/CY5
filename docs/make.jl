using Documenter, CY5
using DocumenterCitations

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
        "Pipeline"   => "pipeline.md",
        "Spaces" => "spaces.md",
        "Example"   => "example.md",
        "Conjectures" => [
            "Overview of Conjectures" => "conj/overview.md",
            "Conjecture 4.2" => "conj/Conj_4_2.md",
            "Conjecture 6.1" => "conj/Conj_6_1.md",
            "Conjecture 7.1" => "conj/Conj_7_1.md",
            "Conjecture 7.2" => "conj/Conj_7_2.md",
            "Conjecture 8.1" => "conj/Conj_8_1.md"
        ],
        "References" => "references.md"
    ],
)
