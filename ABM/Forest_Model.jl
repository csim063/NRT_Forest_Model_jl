"""
This is the primary script of the NRT forest model. This script sets up and runs the model. It can
also record and export data from the model run, and visualise the model. The model can be visualised
as a static image, an interactive image, an image with explanatory plots, or as a video.
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

#* CUSTOM MODULES
include("Modules/Initialise.jl")
include("Modules/Step.jl")

#//-----------------------------------------------------------------------------------------------#
#% DEFINE RUN BEHAVIOUR
n_steps = 100
record_data = false
record_every_n_ticks = 5
export_data = false
export_file_name = "test"
report_time = true
visualisation_type = "static"
video_export_name = "NRT_Video"

#//-----------------------------------------------------------------------------------------------#
#% IMPORT SITE AND DEMOGRAPHIC DATAFRAMES
demography_df = DataFrame(CSV.File("Data/demography.txt"));
site_df = DataFrame(CSV.File("Data/forest.txt"));

#//-----------------------------------------------------------------------------------------------#
#% SETUP MODEL
if record_data == true
        adata = [:pos, :species_ID, :growth_form, :height, :dbh, :age];
end

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
                        planting_frequency = 10,
                        phytothera = false, # Whether to include phytothera disease in the model
                        phyto_global_infection_prob = 0.0001,
                        phyto_local_infection_prob = 0.001,
                        phyto_infectious_radius = 1,
                        phyto_symptoms_dev_prob = 0.1,
                        phyto_mortality_prob = 0.1,
                        );

#//-----------------------------------------------------------------------------------------------#
#% RUN MODEL
if record_data == true
    if report_time == true
        data, _ = @time run!(model, go.agent_step!, go.model_step!, n_steps; 
                        adata = adata, when = [record_every_n_ticks]);
    else
        data, _ = run!(model, go.agent_step!, go.model_step!, n_steps; 
                        adata = adata, when = [record_every_n_ticks]);
    end
    if export_data == true
        CSV.write("Outputs/$export_file_name.csv", data)
    end
else
    if report_time == true
        @time run!(model, go.agent_step!, go.model_step!, n_steps);
    else
        run!(model, go.agent_step!, go.model_step!, n_steps);
    end
end;


#//-----------------------------------------------------------------------------------------------#
#% VISUALISE
#% SETUP PLOTTING VARIABLES-----------------------------------------#
cols = palette("Spectral", nrow(site_df));
speciescolor(a) = a.species_ID == 0 ? :white : cols[a.species_ID];
speciesshape(a) = a.growth_form == 1 ? :circle : :diamond;
speciessize(a) = a.age == 0 ? 8 : (8 + (a.age * 0.025));

plotkwargs = (;
        ac = speciescolor,
        am = speciesshape,
        as = speciessize,
        scatterkwargs = (strokewidth = 0.5,)
);

#% SETUP INTERACTIVE PARAMETER RANGES-------------------------------#
params = Dict(
        :disturbance_freq => 0.1:0.05:1.0,
        :max_disturb_size => 0.1:0.05:1.0,
        :comp_multiplier => 0.0:0.05:3.0,
        :edge_effects => 0:1:1,
        :external_rain => 0:1:1,
        :herbivory => 0:1:1,
        :saplings_eaten => 0:1:1,
        :macro_litter_effect => 0.0:0.01:1.0,
        :ddm => 0:1:1,
        :restoration_planting => 0:1:1,
        :planting_frequency => 1:1:10
);

#% STATIC IMAGE-----------------------------------------------------#
if visualisation_type == "static"
        fig, ax, abmobs  = abmplot(model;
                ac = speciescolor,
                am = speciesshape,
                as = speciessize,
                scatterkwargs = (strokewidth = 0.5,)
        )
        fig
end

#% INTERACTIVE VARIABLES IMAGE--------------------------------------#
if visualisation_type == "interactive"
        using GLMakie
        fig, ax, abmobs = abmplot(model;
                        agent_step! = go.agent_step!,
                        model_step! = go.model_step!,
                        params,
                        plotkwargs...);
        fig
end

#% EXPLORATION PLOT-------------------------------------------------#
#* This plot type includes both the world image but also plots
#* key measures
if visualisation_type == "exploration"
        using Statistics: mean

        ## Agent data to plot
        Tawa(a) = a.species_ID == :1
        Pigeonwood(a) = a.species_ID == :2
        adata = [(Tawa, count), (Pigeonwood, count)]

        ## Global model data to plot
        seedling_count(model) = sum(sum(model.seedlings))
        sapling_count(model) = sum(sum(model.saplings))
        mdata = [seedling_count, sapling_count]

        fig, abmobs = abmexploration(model;
                        agent_step! = go.agent_step!,
                        model_step! = go.model_step!,
                        params,
                        plotkwargs...,
                        adata, alabels =  ["Tawa", "Pigeonwood"],
                        mdata, mlabels = ["Seedlings", "Saplings"]);
        fig
end

#% RECORD A VIDEO---------------------------------------------------#
if visualisation_type == "video"
        abmvideo(
                "Outputs/Videos/$video_export_name.mp4",
                model, go.agent_step!, go.model_step!;
                title ="", frames = 100, framerate = 5,
                showstep = false,
                plotkwargs...
        )
end
