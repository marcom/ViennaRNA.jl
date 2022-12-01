using Test

# show which testset is currently running
showtestset() = println(" "^(2 * Test.get_testset_depth()), "testing ",
                        Test.get_testset().description)

@testset verbose=true "ViennaRNA" begin
    showtestset()
    include("wrappers.jl")
    include("extras.jl")
end
