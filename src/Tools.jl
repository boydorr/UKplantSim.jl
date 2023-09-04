using BritishNationalGrid
using EcoSISTEM.ClimatePref
using AxisArrays
using Unitful
using Distributions
using LinearAlgebra
using EcoSISTEM
using DataFrames
using DataFrameMacros

"""
    OSGR_eastnorth(osgridref::String)

Function to convert between OS grid references and easting/northing
"""
function OSGR_eastnorth(osgridref::String)
    squares = BritishNationalGrid.square_names()
    startingref = findall(osgridref[1:2] .== squares)
    northing = (startingref[1][1] - 1) * 100_000 + parse(Int, osgridref[3:4]) * 1_000
    easting = (startingref[1][2] - 1) * 100_000 + parse(Int, osgridref[5:6]) * 1_000
    [easting, northing]
end

"""
    extractvalues(x::Vector{L},y::Vector{L}, ref::Reference) where L <: Unitful.Length

Function to extract values at specific x,y locations from a `Reference` type of grid.
"""
function extractvalues(x::Vector{L},y::Vector{L}, ref::Reference) where L <: Unitful.Length
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

function extractvalues(x::L, y::L, ref::Reference) where L <: Unitful.Length
    (x <= 7e5m) && (x >= 0.0m) ||
    error("X coordinate is out of bounds")
    (y <= 1.3e6m) && (y >= 0.0m) ||
    error("Y coordinate is out of bounds")
    thisstep = step(ref.array.axes[1].val)
    return ref.array[(x - thisstep/2)..(x + thisstep/2),
              (y - thisstep/2)..(y + thisstep/2)][1,1]
end

"""
    createRef(gridsize::Unitful.Length{Float64}, xmin, xmax, ymin, ymax)

Function to create a reference grid array of type `Reference` with specific x and y bounds.
"""
function createRef(gridsize::Unitful.Length{Float64}, xmin, xmax, ymin, ymax)
    x = length(xmin:gridsize:xmax)
    y = length(ymin:gridsize:ymax)
    refarray = AxisArray(Array{Int64, 2}(undef, Int(floor(x)), Int(floor(y))),Axis{:easting}(xmin:gridsize:xmax),Axis{:northing}(ymin:gridsize:ymax))
    refarray[1:length(refarray)]= collect(1:length(refarray))
    ref = Reference(refarray)
end

LC2015cats = Dict(1.0 => "Broadleaved woodland", 2.0 => "Coniferous woodland", 3.0 => "Arable and horticulture", 4.0 => "Improved grassland", 5.0 => "Neutral grassland", 6.0 => "Calcareous grassland", 7.0 => "Acid grassland", 8.0 => "Fen, marsh and swamp", 9.0 => "Heather", 10.0 => "Heather grassland", 11.0 => "Bog", 12.0 => "Inland rock", 13.0 => "Saltwater", 14.0 => "Freshwater", 15.0 => "Supra-littoral rock", 16.0 => "Supra-littoral sediment", 17.0 => "Littoral rock", 18.0 => "Littoral sediment", 19.0 => "Saltmarsh", 20.0 => "Urban", 21.0 => "Suburban", NaN => "Unknown")
CC2017cats = Dict(1.0 => "Beet", 2.0 => "Field beans", 3.0 => "Grass", 4.0 => "Maize", 5.0 => "Oilseed rape", 6.0 => "Other crops", 7.0 => "Potatoes", 8.0 => "Spring barley", 9.0 => "Spring wheat", 10.0 => "Winter barley", 11.0 => "Winter wheat")

LC2015cats = Dict(1.0 => "Broadleaved woodland", 2.0 => "Coniferous woodland", 3.0 => "Arable and horticulture", 4.0 => "Improved grassland", 5.0 => "Neutral grassland", 6.0 => "Calcareous grassland", 7.0 => "Acid grassland", 8.0 => "Fen, marsh and swamp", 9.0 => "Heather", 10.0 => "Heather grassland", 11.0 => "Bog", 12.0 => "Inland rock", 13.0 => "Saltwater", 14.0 => "Freshwater", 15.0 => "Supra-littoral rock", 16.0 => "Supra-littoral sediment", 17.0 => "Littoral rock", 18.0 => "Littoral sediment", 19.0 => "Saltmarsh", 20.0 => "Urban", 21.0 => "Suburban", NaN => "Unknown", 22.0 => "Beet", 23.0 => "Field beans", 24.0 => "Grass", 25.0 => "Maize", 26.0 => "Oilseed rape", 27.0 => "Other crops", 28.0 => "Potatoes", 29.0 => "Spring barley", 30.0 => "Spring wheat", 31.0 => "Winter barley", 32.0 => "Winter wheat")

"""
    startingArray(bsbi::DataFrame, numspecies::Int64, sf::Int64)

Function to create a grid of starting abundances for ecosystem simulation.
This takes in BSBI locations for a specified number of species `numspecies `
and coarsening factor `sf` for the number of grid squares the abundances should
be spread across.
"""
function startingArray(bsbi::DataFrame, numspecies::Int64, sf::Int64, res::Unitful.Length{Float64} = 1000.0m, xmin::Unitful.Length{Float64} = 500.0m, xmax::Unitful.Length{Float64} = 7e5m, ymin::Unitful.Length{Float64} = 500.0m, ymax::Unitful.Length{Float64} = 1.25e6m)
    ref = createRef(res, xmin, xmax, ymin, ymax)
    fillarray = Array{Int64, 2}(undef, numspecies, length(ref.array))
    ids = sort(unique(bsbi.SppID))
    dict = Dict(zip(ids, 1:length(ids)))
    clustarray = zeros(size(ref.array))
    probarray = zeros(size(ref.array))
    for i in ids
        find_probs!(bsbi, ref, probarray, clustarray, sf, i)
        multi = rand(Multinomial(Int(sum(probarray .> 0) * 1e3), probarray[1:end]))
        fillarray[dict[i], :] .= multi[1:end]
        clustarray .= 0; probarray .= 0
    end
    return fillarray
end

"""
    find_probs!(grouped_tab::JuliaDB.IndexedTable, ref::Reference, probarray::Matrix{Float64}, clustarray::Matrix{Float64}, sf::Int64, spp::Int64)

Function to find probability of abundance across
"""
function find_probs!(bsbi::DataFrame, ref::Reference, probarray::Matrix{Float64}, clustarray::Matrix{Float64}, sf::Int64, spp::Int64)
    refs = filter(g-> g.SppID == spp, bsbi).refval
    xys = convert_coords.(refs, size(ref.array,1))
    newrefs = map(xys) do xy
        newxs = collect(xy[1]:(xy[1] + sf -1)); newys =  collect(xy[2]:(xy[2] + sf-1))
        newxs = newxs[newxs .< 700]; newys = newys[newys .< 1250]
        newrefs = ref.array[newxs, newys][1:end]
        return newrefs
    end
    clustarray[vcat(newrefs...)] .= 1
    identify_clusters!(clustarray)
    for i in unique(clustarray[vcat(newrefs...)])
        clust = findall(clustarray .== i)
        minX =  minimum(map(x-> x[1], clust))
        maxX =  maximum(map(x-> x[1], clust)) + 1
        minY =  minimum(map(y-> y[2], clust))
        maxY =  maximum(map(y-> y[2], clust)) + 1
        H = [maxX-minX 1.0; 1.0 maxY-minY]
        if !isposdef(H) H = [maxX-minX 0.5; 0.5 maxY-minY] end
        mvnorm = MvNormal([mean([maxX, minX]), mean([maxY, minY])], H)
        probarray[clust] .= (map(x -> pdf(mvnorm, [x[1], x[2]]), clust) .* length(clust))
    end
    probarray ./= sum(probarray)
    probarray[probarray .> 0.1] .= 0.1
    probarray ./= sum(probarray)
end

# Function to create clusters from percolated grid
function identify_clusters!(M::AbstractMatrix)
    dimension=size(M)
    # Begin cluster count
    active = findall(M .> 0)
    # Loop through each grid square in M
    Threads.@threads for i in active
        # Find neighbours of M at this location
        neighbours=get_neighbours(M, i[1], i[2], 8)
        inds = convert_coords.(neighbours[:,1], neighbours[:,2], size(M,1))
        already = M[inds] .> 1
        if any(already)
            M[i] = M[neighbours[already,1][1], neighbours[already,2][1]]
        else
            M[i] = maximum(M) + 1
        end
    end
end