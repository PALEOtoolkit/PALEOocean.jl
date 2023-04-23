using Documenter

import PALEOboxes as PB
import PALEOocean

using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src/paleo_references.bib"))

# Collate all markdown files and images folders from PALEOocean/examples/ into a tree of pages
io = IOBuffer()
println(io, "Collating all markdown files from PALEOocean/examples folder:")
examples_folder = "collated_examples"  
examples_path = normpath(@__DIR__, "src", examples_folder)  # temporary folder to collate files
rm(examples_path, force=true, recursive=true)
examples_pages, examples_includes = PB.collate_markdown(
    io, normpath(@__DIR__, "../examples"), @__DIR__, examples_folder;
)
@info String(take!(io))

# include files that load modules etc from PALEOexamples folders
include.(examples_includes)

makedocs(bib, sitename="PALEOocean Documentation", 
    pages = [
        "index.md",
        "Examples and Tutorials" => examples_pages,
        "Design" => [
             "PALEOocean_Domains.md",
        ],
        # no HOWTO docs yes
        "Reference" => [
            "PALEOocean_Reactions.md",
            "PALEOocean_functions.md",
        ],
        "References.md",
        "indexpage.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
)

@info "Local html documentation is available at $(joinpath(@__DIR__, "build/index.html"))"

deploydocs(
    repo = "github.com/PALEOtoolkit/PALEOocean.jl.git",
)
