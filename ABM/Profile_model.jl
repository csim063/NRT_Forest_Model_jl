"""
Script runs the full NRT forest model at a small scale and provides a profiler output for the
`run!()` function (the function that runs the agent and model stepping functions). 

For more details of how profiling works in Julia see 
https://docs.julialang.org/en/v1/manual/profile/
"""

using Agents
using Random

using DelimitedFiles
using DataFrames
using CSV
using Distributions

#* CUSTOM MODULES
include("Modules/Initialise.jl")
include("Modules/Step.jl")

#//------------------------------------------------------------------------------------------------#
#% IMPORT SITE AND DEMOGRAPHIC DATAFRAMES
@info "Importing site and demography data"
demography_df = DataFrame(CSV.File("Data/demography.txt"));
site_df = DataFrame(CSV.File("Data/forest.txt"));

#//------------------------------------------------------------------------------------------------#
#% SETUP AND RUN MODEL INITIALLY TO ALLOW COMPILATION
#*Compilation time only occurs when functions are initially defined and can alter the actual
#* timing and running of the functions thus impacting accurate profiling. To avoid inaccuracy 
#* in the profile results the model is set up and run once with no results captured, before actual
#* profiling occurs

@info "Setting up initial model"
model = Setup.forest_model(forest_area = 4,
                           cell_grain = 4, 
                           n_species = 8,
                           edge_strength = 0.0,
                           max_shade_distance = 32,
                           site_df = site_df,
                           demography_df = demography_df,
                           disturb_freq = 0.000,
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
                           );

@info "Running initial model"
run!(model, go.agent_step!, go.model_step!, 5);

#//------------------------------------------------------------------------------------------------#
#% SETUP MODEL AND PROFILE RUN FUNCTION

@info "Setting up final model"
model = Setup.forest_model(forest_area = 4,
                           cell_grain = 4, 
                           n_species = 8,
                           edge_strength = 0.0,
                           max_shade_distance = 32,
                           site_df = site_df,
                           demography_df = demography_df,
                           disturb_freq = 0.000,
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
                           );

if isinteractive() == false
    @info "Running via console so loading ProfileView"
    using ProfileView
end

@info "Running profile"
@profview run!(model, go.agent_step!, go.model_step!, 1000; showprogress=true)