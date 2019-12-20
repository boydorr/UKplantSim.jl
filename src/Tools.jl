using BritishNationalGrid
using ClimatePref
using AxisArrays
using Unitful
using JuliaDBMeta
using JuliaDB
using Distributions

function OSGR_eastnorth(osgridref::String)
    squares = BritishNationalGrid.square_names()
    startingref = findall(osgridref[1:2] .== squares)
    northing = (startingref[1][1] - 1) * 100_000 + parse(Int, osgridref[3:4]) * 1_000
    easting = (startingref[1][2] - 1) * 100_000 + parse(Int, osgridref[5:6]) * 1_000
    [easting, northing]
end

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

function startingArray(uk::JuliaDB.DIndexedTable, numspecies::Int64)
    ref = createRef(1000.0m, 500.0m, 7e5m, 500.0m, 1.25e6m)
    fillarray = Array{Int64, 2}(undef, numspecies, length(ref.array))
    grouped_tab = @groupby uk (:SppID, :refval) {count = length(:refid)}
    ids = sort(unique(collect(select(uk, :SppID))))
    dict = Dict(zip(ids, 1:length(ids)))
    sppnames = [dict[x] for x in collect(select(grouped_tab, :SppID))]
    refs = collect(select(grouped_tab, :refval))
    counts = collect(select(grouped_tab, :count))
    map(1:length(counts)) do i
        fillarray[sppnames[i], refs[i]] = counts[i] .* 1e3
    end
    return fillarray
end

function startingArray(bsbi::JuliaDB.DIndexedTable, numspecies::Int64, sf::Int64)
    ref = createRef(1000.0m, 500.0m, 7e5m, 500.0m, 1.25e6m)
    fillarray = Array{Int64, 2}(undef, numspecies, length(ref.array))
    grouped_tab = @groupby bsbi (:SppID, :refval) {count = length(:refid)}
    ids = sort(unique(collect(select(bsbi, :SppID))))
    dict = Dict(zip(ids, 1:length(ids)))
    sppnames = [dict[x] for x in collect(select(grouped_tab, :SppID))]
    refs = collect(select(grouped_tab, :refval))
    counts = collect(select(grouped_tab, :count))
    map(1:length(counts)) do i
        x, y = convert_coords(refs[i], size(ref.array, 1))
        xs = collect(x:(x + sf -1)); ys =  collect(y:(y + sf-1))
        xs = xs[xs .< 700]; ys = ys[ys .< 1250]
        newrefs = ref.array[xs, ys][1:end]
        fillarray[sppnames[i], newrefs] .= rand(Multinomial(Int64(counts[i] .* 1e3), length(newrefs)))
    end
    return fillarray
end
function count_neighbours(j::Int64, refs::Vector{Int64}, ref::Reference)
    x,y = Simulation.convert_coords(j, size(ref.array, 1))
    neighbours = Simulation.get_neighbours(Array(ref.array), x, y, 8)
    inds = map((x,y) -> ref.array[x, y], neighbours[:,1], neighbours[:,2])
    return length(inds ∩ newrefs)
end

function startingArray(bsbi::JuliaDB.DIndexedTable, numspecies::Int64, sf::Int64)
    ref = createRef(1000.0m, 500.0m, 7e5m, 500.0m, 1.25e6m)
    fillarray = Array{Int64, 2}(undef, numspecies, length(ref.array))
    grouped_tab = @groupby bsbi (:SppID, :refval) {count = length(:refid)}
    ids = sort(unique(collect(select(bsbi, :SppID))))
    dict = Dict(zip(ids, 1:length(ids)))
    sppnames = [dict[x] for x in collect(select(grouped_tab, :SppID))]
    refs = collect(select(grouped_tab, :refval))
    counts = collect(select(grouped_tab, :count))
    map(1:length(counts)) do i
        x, y = convert_coords(refs[i], size(ref.array, 1))
        xs = collect(x:(x + sf -1)); ys =  collect(y:(y + sf-1))
        xs = xs[xs .< 700]; ys = ys[ys .< 1250]
        newrefs = ref.array[xs, ys][1:end]
        prob = map(j -> count_neighbours(j, Array(newrefs), ref), newrefs)
        if sum(prob) == 0 prob .= 1/length(prob) end
        fillarray[sppnames[i], newrefs] .= rand(Multinomial(Int64(counts[i] .* 1e3), prob./sum(prob)))
    end
    return fillarray
end


refs = select(filter(g-> g.SppID == ids[1], grouped_tab), :refval)
xs,ys = convert_coords(refs, size(ref.array,1))
a = zeros(size(ref.array))
newrefs = map(xs, ys) do x, y
    newxs = collect(x:(x + sf -1)); newys =  collect(y:(y + sf-1))
    newxs = newxs[newxs .< 700]; newys = newys[newys .< 1250]
    newrefs = ref.array[newxs, newys][1:end]
    return newrefs
end
a[vcat(newrefs...)] .= 1
heatmap(a)

identify_clusters!(a)
heatmap(a)

for i in unique(a[vcat(newrefs...)])
    findall(a .== i)
    
end

# Function to create clusters from percolated grid
function identify_clusters!(M::AbstractMatrix)
  dimension=size(M)
  # Begin cluster count
  count=1
  # Loop through each grid square in M
  for x in 1:dimension[1]
    for y in 1:dimension[2]

      # If square is marked as 1, then apply cluster finding algorithm
      if M[x,y]==1.0
        # Find neighbours of M at this location
        neighbours=get_neighbours(M, x, y, 8)
        # Find out if any of the neighbours also have a value of 1, thus, have
        # not been assigned a cluster yet
        cluster = vcat(mapslices(x->M[x[1],x[2]] .== 1, neighbours, dims=2)...)
        # Find out if any of the neighbours have a value > 1, thus, have already
        # been assigned a cluster
        already=vcat(mapslices(x->M[x[1],x[2]] .> 1, neighbours, dims=2)...)
        # If any already assigned neighbours, then assign the grid square to this
        # same type
          if any(already)
            neighbours=neighbours[already,:]
            M[x,y]=M[neighbours[1,1],neighbours[1,2]]
          # If none are assigned yet, then create a new cluster
          else
            count=count+1
            neighbours=neighbours[cluster,:]
            M[x,y]=count
            map(i->M[neighbours[i,1],neighbours[i,2]]=count, 1:size(neighbours,1))
        end
      end
    end
  end
end
