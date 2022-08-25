"""
## TODO Finish writing overall script description
This script acts as the main or consolidation script for the model. In this
script all the required functions from ./Modules/ are called
"""

using Agents
using Random

using DelimitedFiles
using DataFrames
using CSV
using Distributions

using Colors
using ColorBrewer

using InteractiveDynamics
using CairoMakie

## Custom modules
include("Modules/Initialise.jl")

#//---------------------------------------------------------------------------#
#% Import data from files
demography_df = DataFrame(CSV.File("Data/demography.txt"))
site_df = DataFrame(CSV.File("Data/forest.txt"))

#//---------------------------------------------------------------------------#
#% Setup model
model = Setup.forest_model(forest_area = 4,
                           cell_grain = 4, 
                           edge_strength = 1.0,
                           max_shade_distance = 32,
                           site_df = site_df,
                           demography_df = demography_df
                           )

#//---------------------------------------------------------------------------#
#% ...


#//---------------------------------------------------------------------------#
#% Visualise
cols = palette("Spectral", nrow(site_df));
speciescolor(a) = a.species_ID == 0 ? :white : cols[a.species_ID]

speciesshape(a) = a.growth_form == 1 ? :circle : :diamond
speciessize(a) = a.age == 0 ? 8 : (8 + (a.age * 0.025))

fig,df  = abmplot(model;
        ac = speciescolor,
        am = speciesshape,
        as = speciessize
)
fig

#//---------------------------------------------------------------------------#
#% Export data
