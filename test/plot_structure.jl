import Luxor

@testset "plot_structure" begin
    w = "((((.....))))."
    s = "GGCGAAUACCGCCU"
    @test plot_structure(w) isa Luxor.Drawing
    @test plot_structure(w; sequence=s) isa Luxor.Drawing
    for layout in [:simple, :naview, :circular, :turtle, :puzzler]
        @test plot_structure(w; sequence=s, layout_type=layout) isa Luxor.Drawing
        @test plot_structure(w; sequence=s, layout_type=layout,
                             base_colors=rand(length(w))) isa Luxor.Drawing        
    end
end
