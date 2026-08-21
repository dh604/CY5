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
    format   = Documenter.HTML(edit_link = "main", mathengine = MathJax3()),
    pages    = [
        "Home"       => "index.md",
        "Gromov-Witten computations" => [
            "The package `GW_CY5.jl`" => "cy5_overview.md",
            "Pipeline"   => "pipeline.md",
            "Spaces" => "spaces.md",
            "Examples"   => [
                    "Julia code in the paper" => "code_from_paper.md",
                    "Details of the Pipeline" => "example.md"
                ],
            "Evidence" => [
                "Overview of GW evidence" => "conj/overview.md",
                " \$X_{\\text{3d-strip}}\\times\\mathbb{C}^2\$ (Evidence 2.8)" => "conj/3d_strip_xC2.md",
                "Closed vertex (Evidence 2.13)" => "conj/closed_vertex.md",
                "Gauge theory (Evidence 3.5)" => "conj/gw_gauge.md",
                "Local \$\\mathbb{P}^2\$ (Evidence 4.2)" => "conj/P2_111.md",
                "Local \$\\mathbb{P}^3\$ (Evidence 4.8)" => "conj/P3_a_4a.md",
                "Local \$\\mathbb{P}^2\$, \$\\mathbb{P}^1\\times\\mathbb{P}^1\$ (Evidence 5.5)" => "conj/cy5_from_cy4.md"
            ]
        ],
        "Pandharipande-Thomas computations" => [
            "Overview" => "pt/pt_cy5_overview.md",
            "Example" => "pt/pt_example.md",
            "Evidence" => [
                "Overview of PT evidence" => "pt/pt_conj/overview.md",
                "Closed vertex (Evidence 2.12)" => "pt/pt_conj/closed_vertex.md",
                "Local \$\\mathbb{P}^2\$ (Evidence 4.4 & 5.7)" => "pt/pt_conj/local_P2.md",
                "Local \$\\mathbb{P}^3\$ (Evidence 4.9)" => "pt/pt_conj/local_P3.md",
                "Local \$\\mathbb{P}^1 \\times \\mathbb{P}^1\$ (Evidence 5.7)" => "pt/pt_conj/local_P1xP1.md"
            ]
        ],
        "References" => "references.md"
    ];
    doctest=true
)
