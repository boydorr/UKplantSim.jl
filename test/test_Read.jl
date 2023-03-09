using Test
using UKplantSim

file = "../data/CEH_landcover_2015.tif"
@test_nowarn readLC(file)
file = "../data/Crop2017.tif"
@test_nowarn readCrop(file)
file = "../data/HuttonSoils.tif"
@test_nowarn readSoils(file)