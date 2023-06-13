import Aqua
using ViennaRNA

@testset "Aqua.test_all" begin
    showtestset()
    Aqua.test_all(ViennaRNA)
end
