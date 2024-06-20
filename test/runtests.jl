

using Test
using Documenter

ENV["GKSwstype"] = "100" # to run Plots.jl GR on a headless system

@testset "PALEOocean all" begin

@testset "PALEOocean test" begin
    include("runocean3boxtests.jl")

    # include("runoceanMITgcmtests.jl")
end

@testset "PALEOocean/examples" begin

include("../examples/transport_examples/runtests.jl")

include("../examples/ocean3box/runtests.jl")

include("../examples/PTBClarkson2014/runtests.jl")

include("../examples/romglb/runtests.jl")

include("../examples/blacksea/runtests.jl")

end

doctest(PALEOocean; manual=false)  

end

delete!(ENV, "GKSwstype"); # undo workaround for Plots.jl on a headless system