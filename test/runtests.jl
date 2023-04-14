

using Test
using Documenter


@testset "PALEOocean all" begin

include("../examples/transport_examples/runtests.jl")

doctest(PALEOocean; manual=false)  

end