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

## Custom modules
include("Modules/Initialise.jl")

#//---------------------------------------------------------------------------#
#% Import data from files
demography_df = DataFrame(CSV.File("Data/demography.txt"))
site_df = DataFrame(CSV.File("Data/forest.txt"))


#//---------------------------------------------------------------------------#
#% Setup model
Setup.setup() #*Note this is how to call the functions in the modules in the files you have included

#//---------------------------------------------------------------------------#
#% ...


#//---------------------------------------------------------------------------#
#% Visualise
## TODO fix this
using InteractiveDynamics
using CairoMakie
model = Setup.forest_model(site_df = site_df)

fig, _ = abmplot(model;)
fig

#//---------------------------------------------------------------------------#
#% Export data
