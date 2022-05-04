using Statistics: mean
using ColorSchemes: ColorScheme, colorschemes
using Luxor: @draw, background, circle, fontsize, line, Point, sethue, text

# TODO
# - text font for bases
# - legend for base_colors gradient
# - choose png or svg output
# - be able to ignore base_colors if not desired,
#   with base_colorscheme == flag_id it is white by default
# - make `constants` block configurable

"""
    plot_structure(structure; [sequence, layout_type, base_colors, base_colorscheme])

Plot a secondary structure to a PNG image.
"""
function plot_structure(structure::AbstractString;
                        sequence::AbstractString=" "^length(structure),
                        layout_type::Symbol=:simple,
                        base_colors::Vector=zeros(length(structure)),
                        base_colorscheme::ColorScheme=colorschemes[:flag_id])
    # constants
    base_radius = 10
    font_size = 20
    background_color = "white"
    base_circle_color = "black"
    base_text_color = "black"
    backbone_color = "black"
    basepair_color = "blue"
    scale_plot_coords = if layout_type == :circular
        200.0
    elseif layout_type == :turtle || layout_type == :puzzler
        1.3
    else
        2.5
    end
    x_pad = 5
    y_pad = 5

    function line_connect_bases(xa, ya, xb, yb, color)
        sethue(color)
	# TODO: line goes over the edge of a base, this is fixed in
	# draw_base by first overwriting with background color
	line(Point(xa, ya), Point(xb, yb), :stroke)
    end
    draw_backbone(xa, ya, xb, yb) =
	line_connect_bases(xa, ya, xb, yb, backbone_color)
    draw_basepair(xa, ya, xb, yb) =
	line_connect_bases(xa, ya, xb, yb, basepair_color)
    function draw_base(x, y, txt::AbstractString, fillcolor)
        p = Point(x, y)
        # overwrite with background
	sethue(background_color)
	circle(Point(x, y), base_radius, :fill)
        # base_color
        sethue(fillcolor)
	circle(p, base_radius, :fill)
        # base circle
	sethue(base_circle_color)
	circle(p, base_radius, :stroke)
        # text
        sethue(base_text_color)
	text(txt, p, halign=:center, valign=:middle)
    end

    n = length(structure)
    # TODO: we need the structure as a String for plot_coords,
    # and as a Pairtable to find basepairs
    pt = Pairtable(structure)
    length(sequence) == n ||
        throw(ArgumentError("structure and sequence must have same length"))
    length(base_colors) == n ||
        throw(ArgumentError("base_colors vector must have same length as structure"))

    xs, ys = plot_coords(structure; plot_type=layout_type)
    # center coords, and calculate width and height of image
    # TODO: not sure if optimal this way, seems ok-ish
    xs .*= scale_plot_coords
    ys .*= scale_plot_coords
    x_origin = mean(xs)
    y_origin = mean(ys)
    xs .-= x_origin
    ys .-= y_origin
    x_min, x_max = extrema(xs)
    y_min, y_max = extrema(ys)
    x_width = 2 * (ceil(max(abs(x_min), abs(x_max))) + base_radius + x_pad)
    y_width = 2 * (ceil(max(abs(y_min), abs(y_max))) + base_radius + y_pad)

    # Note:
    # in Luxor, (0,0) is the center,
    #           (-width/2, -height/2) is top left
    #       and (width/2, height/2) is lower right
    @draw begin
	background(background_color)
	fontsize(font_size)
	# backbone
	for i = 1:n
	    if i != n
                draw_backbone(xs[i], ys[i], xs[i+1], ys[i+1])
	    end
	end
	# basepairs
	for i = 1:n
	    j = pt[i]
	    (j == 0 || i > j) && continue
            draw_basepair(xs[i], ys[i], xs[j], ys[j])
	end
        # bases (these must come last, as we overwrite too long
        # backbone and basepair lines here)
	for i = 1:n
            draw_base(xs[i], ys[i], string(sequence[i]),
                      base_colorscheme[base_colors[i]])
	end
    end x_width y_width
end

plot_structure(pt::Pairtable;
               sequence::AbstractString=" "^length(pt),
               layout_type::Symbol=:naview,
               base_colors::Vector=zeros(length(pt)),
               base_colorscheme::ColorScheme=colorschemes[:flag_id]) =
                   plot_structure(String(pt);
                                  sequence, base_colors, base_colorscheme)


using CairoMakie


"""
    plot_structure(structure; [sequence, targetdir, layout_type, colorscheme])

Plot a secondary structure to a PNG image or PDF file depending on targetdir ending.
"""
function plot_structure_makie(
    structure::AbstractString;
    sequence::AbstractString=" "^length(structure),
    filepath::String = "",
    layout_type::Symbol=:simple,
    colorscheme::Symbol=:lightrainbow)

    fc = FoldCompound(sequence)
    partfn(fc)
    pt = Pairtable(structure)
    x_coords, y_coords = plot_coords(structure; plot_type=layout_type)

    markersize = 100 / sqrt(length(structure))
    positions = [(a, b) for (a, b) in zip(x_coords, y_coords)]
    pairs = [(i, pt[i]) for i = 1:length(pt) if i < pt[i]]
    f = Figure()
    ax = Axis(f[1, 1])
    xlims!(
        round(Int, minimum(x_coords)) - 20,
        round(Int, maximum(x_coords)) + 20)
    ylims!(
        round(Int, minimum(y_coords)) - 20,
        round(Int, maximum(y_coords)) + 20)
    hidedecorations!(ax)
    hidespines!(ax)

    for (i, j) in pairs
        lines!(
            [x_coords[i], x_coords[j]],
            [y_coords[i], y_coords[j]],
            color = :black,
            linestyle = :dot,
            linewidth = 1,
        )
    end
    scatterlines!(ax, x_coords, y_coords,
        markersize = markersize,
        markercolor = prob_of_basepairs(fc, pt),
        markercolormap = colorscheme,
        markercolorrange = (0, 1), # probabilities [0, 1]
        linewidth = 3,
    )
    text!(string.(collect(sequence)),
        position = positions,
        align = (:center, :center),
        textsize = markersize / 2,
    )
    Colorbar(
        f[2, 1],
        vertical = false,
        colormap = colorscheme,
        width = 500,
        height = 12,
    )
    ax.aspect = DataAspect()
    if ! isempty(filepath)
        save(filepath, f)
    end
    return f
end
