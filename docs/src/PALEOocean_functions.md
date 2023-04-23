# PALEOocean functions

Helper functions for use by Reactions.

```@meta
CurrentModule = PALEOocean.Ocean
```
## Configuring Domains and Variables

```@docs
set_model_domains
find_transport_vars
```

## Constructing and using transport matrices
```@docs
add_loop!
prepare_transport
do_transport
do_transport_tr
```

### Optimized transport using transport matrices with a common sparsity pattern
```@docs
TrsptCSC
create_common_sparsity_tr!
```

### Optimized transport using SIMD packed vectors
```@docs
PackedBuffer
```
