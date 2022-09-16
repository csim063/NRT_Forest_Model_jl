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
include("Modules/Step.jl")

#//---------------------------------------------------------------------------#
#% Import data from files
demography_df = DataFrame(CSV.File("Data/demography.txt"))
site_df = DataFrame(CSV.File("Data/forest.txt"))

#//---------------------------------------------------------------------------#
#% Setup model
## First setup the dataframe where data will be saved
adata = [:pos, :species_ID, :growth_form, :height, :dbh, :age]

model = Setup.forest_model(forest_area = 4,
                           cell_grain = 4, 
                           n_species = 8,
                           edge_strength = 1.0,
                           max_shade_distance = 32,
                           site_df = site_df,
                           demography_df = demography_df,
                           disturb_freq = 0.200,
                           max_disturb_size = 0.40,
                           comp_multiplier = 1.60,
                           edge_effects = false,
                           external_rain = false,
                           ext_dispersal_scenario = "equal",
                           herbivory = false,
                           saplings_eaten = false,
                           macro_litter_effect = 0.10,
                           ddm = false,
                           restoration_planting = false,
                           planting_frequency = 10
                           )

#//---------------------------------------------------------------------------#
#% Step model
#step!(model, go.agent_step!, go.model_step!, 1)
data, _ = run!(model, go.agent_step!, go.model_step!, 5; adata = adata, when = [5])
data[1:100, :]

# Quick and dirty performance benchmarking
@time run!(model, go.agent_step!, go.model_step!, 5);

#//---------------------------------------------------------------------------#
#% Visualise
cols = palette("Spectral", nrow(site_df));
speciescolor(a) = a.species_ID == 0 ? :white : cols[a.species_ID]

speciesshape(a) = a.growth_form == 1 ? :circle : :diamond
speciessize(a) = a.age == 0 ? 8 : (8 + (a.age * 0.025))

fig, _  = abmplot(model;
        ac = speciescolor,
        am = speciesshape,
        as = speciessize
)
fig

#//---------------------------------------------------------------------------#
#% Export data
