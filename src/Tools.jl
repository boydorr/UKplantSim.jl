using BritishNationalGrid
using ClimatePref
using AxisArrays
using Unitful

function OSGR_eastnorth(osgridref::String)
    squares = BritishNationalGrid.square_names()
    startingref = findall(osgridref[1:2] .== squares)
    northing = (startingref[1][1] - 1) * 100_000 + parse(Int, osgridref[3:4]) * 1_000
    easting = (startingref[1][2] - 1) * 100_000 + parse(Int, osgridref[5:6]) * 1_000
    [easting, northing]
end

function extractvalues(x::Vector{typeof(1m)},y::Vector{typeof(1m)},
    ref::Reference)
    all(x .<= 7e5m) && all(x .>= 0.0m) ||
    error("X coordinate is out of bounds")
    all(y .<= 1.3e6m) && all(y .>= 0.0m) ||
    error("Y coordinate is out of bounds")
    thisstep = step(ref.array.axes[1].val)
    map(x , y) do lon, lat
        return ref.array[(lon - thisstep/2)..(lon + thisstep/2),
              (lat - thisstep/2)..(lat + thisstep/2)][1,1]
    end
end

function extractvalues(x::typeof(1m),y::typeof(1m),
    ref::Reference)
    (x <= 7e5m) && (x >= 0.0m) ||
    error("X coordinate is out of bounds")
    (y <= 1.3e6m) && (y >= 0.0m) ||
    error("Y coordinate is out of bounds")
    thisstep = step(ref.array.axes[1].val)
    return ref.array[(x - thisstep/2)..(x + thisstep/2),
              (y - thisstep/2)..(y + thisstep/2)][1,1]
end

"""
    create_reference(gridsize::Float64, xmin, xmax, ymin, ymax)

Function to create a reference grid array of type `Reference`.
"""
function createRef(gridsize::Unitful.Length{Float64}, xmin, xmax, ymin, ymax)
    x = length(xmin:gridsize:xmax)
    y = length(ymin:gridsize:ymax)
    refarray = AxisArray(Array{Int64, 2}(undef, Int(floor(x)), Int(floor(y))),Axis{:easting}(xmin:gridsize:xmax),Axis{:northing}(ymin:gridsize:ymax))
    refarray[1:length(refarray)]= collect(1:length(refarray))
    ref = Reference(refarray)
end

LC2015cats = Dict(1.0 => "Broadleaved woodland", 2.0 => "Coniferous woodland", 3.0 => "Arable and horticulture", 4.0 => "Improved grassland", 5.0 => "Neutral grassland", 6.0 => "Calcareous grassland", 7.0 => "Acid grassland", 8.0 => "Fen, marsh and swamp", 9.0 => "Heather", 10.0 => "Heather grassland", 11.0 => "Bog", 12.0 => "Inland rock", 13.0 => "Saltwater", 14.0 => "Freshwater", 15.0 => "Supra-littoral rock", 16.0 => "Supra-littoral sediment", 17.0 => "Littoral rock", 18.0 => "Littoral sediment", 19.0 => "Saltmarsh", 20.0 => "Urban", 21.0 => "Suburban", NaN => "Unknown")
