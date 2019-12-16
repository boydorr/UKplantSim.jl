using NetCDF
using UKclim
using MyUnitful
using Unitful
using Unitful.DefaultSymbols
cd("/home/claireh/Documents/UK")
ncinfo("CHESS/prec/chess_precip_201001.nc")


times = collect(2010year:1month:2011year+9month)
dir = "CHESS/prec1/"
sun = readCHESS(dir, "precip", times)
