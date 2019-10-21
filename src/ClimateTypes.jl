using AxisArrays
using Unitful
using RecipesBase

import ClimatePref.AbstractClimate
import AxisArrays.axes
"""
    HadUK <: AbstractClimate

Type that houses data extracted from HadUK raster files.
"""
mutable struct HadUK <: AbstractClimate
    array::AxisArray
    function HadUK(array::AxisArray)
        typeof(collect(axes(array, 3).val)[1])<: Unitful.Time ||
            error("Third dimension of array must be time")
        new(array)
    end
end
@recipe function f(hd::HadUK, tm::Unitful.Time)
    background_color   := :lightblue
    background_color_outside := :white
    grid := false
    aspect_ratio := 1
    seriestype  :=  :heatmap
    colorbar_title := string(unit(hd.array[1]))
    transpose(ustrip.(hd.array[:, :, tm]))
end

"""
    LandCover <: AbstractClimate

Type that houses data extracted from CEH landcover 2015 rasters.
"""
mutable struct LandCover <: AbstractClimate
    array::AxisArray
    function LandCover(array::AxisArray)
        length(size(array)) == 2 ||
            error("Landcover is two-dimensional")
        new(array)
    end
end

@recipe function f(lc::LandCover)
    color := [:red2,:chartreuse4, :red4, :chartreuse, :lightgreen, :yellowgreen, :gold4, :yellow, :darkmagenta, :lightpink,  :seagreen, :lavender, :navy, :blue, :gold3,:gold3, :khaki1, :khaki1, :lightslateblue, :black, :grey]
    background_color   := :lightblue
    background_color_outside := :white
    grid := false
    aspect_ratio := 1
    seriestype  :=  :heatmap
    transpose(lc.array)
end
