# Skeleton ocean-atmosphere configuration

A "skeleton" ocean-atmosphere configuration, defining Domains and fluxes but no biogeochemistry. Illustrative fluxes are configured for a C, O, P model.

## yaml configuration file
The model configuration (file `examples/skeleton_configuration/PALEO_examples_oceanskeleton_cfg.yaml`) contains:
```@eval
import Markdown
Markdown.parse(
    """```julia
    $(read(joinpath(ENV["PALEO_EXAMPLES"], "skeleton_configuration/PALEO_examples_oceanskeleton_cfg.yaml"), String))
    ```"""
)
```

## Variables defined

    julia> include("PALEO_examples_oceanskeleton.jl")

```@setup oceanskeleton
include(joinpath(ENV["PALEO_EXAMPLES"], "skeleton_configuration/PALEO_examples_oceanskeleton.jl")) # hide
```
```@example oceanskeleton
show(PB.show_variables(model); allrows=true, allcols=true, eltypes=false, show_row_number=false)
```

ocean, oceansurface, and oceanfloor Variables are standard grid variables, provided by the ocean transport Reaction
(a `ReactionOceanTransport3box` in this case).

Flux target Variables to accumulate exchange fluxes are provided by `ReactionFluxTarget` Reactions, the illustrative
configuration here defines fluxes for a C, O, P model with atmospheric state variables CO2, O2, and ocean state
variables including DIC, TAlk, O2 and P.