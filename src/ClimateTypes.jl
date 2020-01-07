using AxisArrays
using Unitful
using RecipesBase

import ClimatePref.AbstractClimate
import AxisArrays.axes
import Plots.cgrad
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
    CHESS <: AbstractClimate

Type that houses data extracted from CHESS raster files.
"""
mutable struct CHESS <: AbstractClimate
    array::AxisArray
    function CHESS(array::AxisArray)
        typeof(collect(axes(array, 3).val)[1])<: Unitful.Time ||
            error("Third dimension of array must be time")
        new(array)
    end
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
    color := cgrad([:red2,:chartreuse4, :red4, :chartreuse, :lightgreen, :yellowgreen, :gold4, :yellow, :darkmagenta, :lightpink,  :seagreen, :lavender, :navy, :blue, :gold3,:gold3, :khaki1, :khaki1, :lightslateblue, :black, :grey])
    background_color   := :lightblue
    background_color_outside := :white
    grid := false
    aspect_ratio := 1
    seriestype  :=  :heatmap
    transpose(lc.array)
end

@recipe function f(lc::Dict)
    order = sortperm(Float64.(keys(lc)))
    lcnames = [LC2015cats[x] for x in Float64.(keys(lc))[order]]
    colorder = Float64.(keys(lc))[order]
    colorder[end] = 22
    color := cgrad([:red2,:chartreuse4, :red4, :chartreuse, :lightgreen, :yellowgreen, :gold4, :yellow, :darkmagenta, :lightpink,  :seagreen, :lavender, :navy, :blue, :gold3,:gold3, :khaki1, :khaki1, :lightslateblue, :black, :grey, :white][Int.(colorder)])
    grid := false
    seriestype  :=  :bar
    legend := false
    size := (1400, 1000)
    xrotation := 55
    guidefontsize := 12
    tickfontsize := 12
    xticks := :all
    lcnames, Int.(values(lc))[order]
end

"""
    CropCover <: AbstractClimate

Type that houses data extracted from CEH crop cover 2017 rasters.
"""
mutable struct CropCover <: AbstractClimate
    array::AxisArray
    function CropCover(array::AxisArray)
        length(size(array)) == 2 ||
            error("Landcover is two-dimensional")
        new(array)
    end
end

@recipe function f(cc::CropCover)
    color := cgrad([:grey, :salmon, :red, :seagreen, :chartreuse, :yellow, :lightblue, :royalblue, :tan, :wheat, :darkorange,  :orange4])
    background_color   := :lightblue
    background_color_outside := :white
    grid := false
    aspect_ratio := 1
    seriestype  :=  :heatmap
    transpose(cc.array)
end
@recipe function f(cc::Dict)
    order = sortperm(Float64.(keys(cc)))
    ccnames = [CC2017cats[x] for x in Float64.(keys(cc))[order]]
    colorder = Float64.(keys(cc))[order]
    colorder[end] = 22
    color := cgrad([:grey, :salmon, :red, :seagreen, :chartreuse, :yellow, :lightblue, :royalblue, :tan, :wheat, :darkorange,  :orange4][Int.(colorder)])
    grid := false
    seriestype  :=  :bar
    legend := false
    size := (1400, 1000)
    xrotation := 55
    guidefontsize := 12
    tickfontsize := 12
    xticks := :all
    ccnames, Int.(values(cc))[order]
end
