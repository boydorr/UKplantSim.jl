using Unitful
using Unitful.DefaultSymbols
using MyUnitful
using AxisArrays
using NetCDF
using Compat

import Unitful.°, Unitful.°C, Unitful.mm
import ArchGDAL
import Base.read
const AG = ArchGDAL

unitdict = Dict("K" => K, "m" => m, "J m**-2" => J/m^2, "m**3 m**-3" => m^3, "degC" => °C)

"""
    searchdir(path,key)

Function to search a directory `path` using a given `key` string.
"""
searchdir(path,key) = filter(x->Compat.occursin(key, x), readdir(path))
dir = "data/HadUK/tas/"
param = "tas"

"""
    read(f, filename)

Function to read raster file into julia.
"""
function read(f, filename)
    return AG.registerdrivers() do
        AG.read(filename) do dataset
            f(dataset)
        end
    end
end

"""
    readlc(file::String)

Function to import a selected CEH land cover file from a path string. Optional arguments for extent if not Great Britain.
"""
function readLC(file::String, GB::Bool=true)
    if GB
        xmin = 0, xmax = 7e5, ymin = 0, ymax = 1.3e6
    else
        xmin = 1.8e5, xmax = 3.7e5, ymin = 3e5, ymax = 4.6e5
    end
    txy = [Float64, Int64(1), Int64(1), Float64(1)]
    #
    read(file) do dataset
        #txy[1] = AG.getdatatype(AG.getband(dataset, 1))
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        txy[4] = AG.getnodatavalue(AG.getband(dataset, 1))
        print(dataset)
    end

    a = Array{txy[1], 2}(Compat.undef, txy[2], txy[3])
    read(file) do dataset
        bd = AG.getband(dataset, 1);
        AG.read!(bd, a);
    end;
    lat, long = size(a, 1), size(a, 2);
    step = abs(xmin - xmax) / lat;
    latitude = (xmin+ step):step:xmax
    longitude = (ymin+ step):step:ymax
    size(longitude)
    lc = AxisArray(a[:, end:-1:1], Axis{:latitude}(latitude), Axis{:longitude}(reverse(longitude)));

    if txy[1] <: AbstractFloat
        lc[lc .== lc[1]] *= NaN;
    end;
    return LandCover(lc)
end

"""
    readHadUK(file::String)

Function to import HadUK data into Julia from particular parameter.
"""
function readHadUK(dir::String, param::String, times::Vector{T}) where T<: Unitful.Time
    files = searchdir(dir, ".nc")
    lat = ncread(joinpath(dir, files[1]), "projection_y_coordinate")
    lon = ncread(joinpath(dir, files[1]), "projection_x_coordinate")
    units = ncgetatt(joinpath(dir, files[1]), param, "units")
    units = unitdict[units]
    array = map(files) do f
        ncread(joinpath(dir, f), param)
    end
    array = cat(dims = 3, array...)
    array[array .≈ ncgetatt(joinpath(dir, files[1]), param, "_FillValue")] .= NaN
    array = array * 1.0 * units

    # If temperature param, need to convert from Kelvin
    if typeof(units) <: Unitful.TemperatureUnits
        array = uconvert.(K, array)
    end
    uk = AxisArray(array, Axis{:longitude}(lon * m), Axis{:latitude}(lat * m), Axis{:month}(times))
    return HadUK(uk[0.0m..1e6m, 0.0m..1.25e6m, :])
end
