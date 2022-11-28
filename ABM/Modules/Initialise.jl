"""
Module primarily contains the function to defiine the agents and setup the model for the ABM. 
Additional miscellaneous helper functions only useful to the setup of the model are also contained
in this module. The code contained in this module can be thought of as equivalent to the setup
function/button of NetLogo.
"""
module Setup
    using Agents
    using Random
    using StatsBase
    using DataFrames

    #* CUSTOM MODULES
    include("Helper_functions.jl")
    include("Demographic_assignments.jl")

    #//-------------------------------------------------------------------------------------------#
    #% DEFINE AGENTS
    
    # This specialised function is defined by Agents.jl itself (see 
    # https://juliadynamics.github.io/Agents.jl/stable/api/#Agents.@agent). The function creates and
    # defines the agents, referred to as Trees, for the model. They are defined as GridAgents{2} 
    # meaning they must be placed and behave on a 2D grid model. Agents have two inbuilt properties 
    # and seven custom properties. These are:
    # - `id::Int`: Unique agent ID.
    # - `pos::NTuple{2, Int}`: Coordinates of agent.
    # - `species_ID::Int`: Value indicating what tree species the agent is.
    # - `patch_here_ID::Int`: Patch ID of the cell the agent is currently on.
    # - `growth_form::Int`: Growth type of agent, 1 = trees; 2 = tree ferns.
    # - `height::Float64`: Height in meters of agent.
    # - `dbh::Float64`: Diameter at breast height in meters of agent.
    # - `age::Float64`: Age (number of ticks) agent has been alive
    # - `previous_growth::Array`: Growth penalty list for the last 5 ticks worth of growth 
    # - `phytothera_infected::Bool`: Whether the agent is infected by soil borne pathogens (e.g. Phytophthora)
    # - `phytothera_infection_age::Int`: Age (number of ticks) since the agent has been infected by soil borne pathogens (e.g. Phytophthora)
    # - `phytothera_symtomatic::Bool`: Whether the agent is showing symptoms of soil borne pathogens (e.g. Phytophthora)

    @agent Tree GridAgent{2} begin 
        species_ID::Int64
        patch_here_ID::Int64
        growth_form::Int64
        height::Float64
        dbh::Float64
        age::Float64
        previous_growth::Array
        phytothera_infected::Bool
        phytothera_infected_age::Int64
        phytothera_symptomatic::Bool
    end

    #//-------------------------------------------------------------------------------------------#
    #% DEFINE WORLD
    """
    This is the primary setup function. It defines all the values of the model, including the 
    dimensions and properties of the space, the inital parameter values for all cells and agents, 
    and even the scheduler for undertaking agent behaviours.

    User inputs:
    - `forest_area::Int64`: Size in hectares of the modelled habitat patch
    - `cell_grain::Int64`: Size in meters of an individual cell in the model (cell_grain x cell_grain)
    - `n_species::Int64`: Number of species modelled
    - `edge_strength::Float64`: Value which determines over what distance edge effects occur (0-1)
    - `max_shade_distance::Int64`: Maximum distance over which shading from neighbours can occur
    - `site_df::DataFrame`: Data defining the initial site properties to be modelled
    - `demography_df::DataFrame`: Species specific demography values for the species to be modelled
    - `seed::Int64`: Random seed value to allow for repeatability
    - `disturb_freq::Float64`: Probaility of a disturbance occurring in a tick
    - `max_disturb_size::Float64`: Maximum proportion of patch that may be impacted by a single disturbance
    - `comp_multiplier::Float64`: Constant impacting the strength of competition that each species 
            experiences. This constant is α in α * (height / shade_height) ^ 0.5 based on 
            Dislich et *al* 2009.
    - `edge_effects::Bool`: Whether to include for edge effects or not
    - `external_rain::Bool`: Whether to have seeds disperse into patch from beyond the patch
    - `ext_dispersal_scenario::String`: Whether to assign external seeds to each species equally
        ("equal") or by abundance ("abundance")
    - `herbivory::Bool`: Whether to include for herbivory effects or not
    - `saplings_eaten::Bool`: Whether saplings are impacted by herbivory or not
    - `macro_litter_effect::Float64`: Probability of a sapling being killed by a macro-litter fall
    - `ddm::Bool`: Whether to include density dependent mortality or not
    - `restoration_planting::Bool`: Whether to include restoration planting or not
    - `planting_frequency::Int64`: How often (how many ticks) does restoration planting occur
    - `phytothera::Bool`: Whether to include phytothera disease or not
    - `phyto_global_infection_prob::Float64`: Probability of a tree being infected by phytothera due
        to global chance
    - `phyto_local_infection_prob::Float64`: Probability of a tree being infected by phytothera due
        to local chance from an infected neighbour
    - `phyto_infectious_radius::Int64`: Radius over which phytothera can be spread from an infected
        tree
    - `phyto_symptom_prob::Float64`: Probability of a tree developing symptoms of phytothera in any
        given tick
    - `phyto_mortality_prob::Float64`: Probability of a tree dying from phytothera in any given tick
    """
    function forest_model(;
        forest_area::Int64 = 16,
        cell_grain::Int64 = 4,
        n_species::Int64 = 8,
        edge_strength::Float64 = 0.0,
        max_shade_distance::Int64 = 32,
        site_df::DataFrame = site_df,
        demography_df::DataFrame = demography_df,
        seed::Int64 = 999,
        disturb_freq::Float64 = 0.100,
        max_disturb_size::Float64 = 0.40,
        comp_multiplier::Float64 = 1.60,
        edge_effects::Bool = false,
        external_rain::Bool = false,
        ext_dispersal_scenario::String = "equal",
        herbivory::Bool = false,
        saplings_eaten::Bool = false,
        macro_litter_effect::Float64 = 0.10,
        ddm::Bool = false,
        restoration_planting::Bool = false,
        planting_frequency::Int64 = 10,
        phytothera::Bool = false,
        phyto_global_infection_prob::Float64 = 0.0001,
        phyto_local_infection_prob::Float64 = 0.001,
        phyto_infectious_radius::Int64 = 1,
        phyto_symptoms_dev_prob::Float64 = 0.1,
        phyto_mortality_prob::Float64 = 0.1,
        )


        ###--------------------------------DEFINE SPACE--------------------------------###
        dims = trunc(Int, (sqrt(forest_area * 1e4)) / cell_grain)

        space = GridSpaceSingle((dims, dims); periodic = false, metric = :chebyshev);
        rng = MersenneTwister(seed)

        ###-------------------------DEFINE PROPERTY VARIABLES-------------------------###
        seedling_survival = demography_df.seedling_survival
        sapling_survival = demography_df.sapling_survival
        seedling_transition = demography_df.seedling_transition
        seedling_mortality = 1 .- (seedling_survival .+ seedling_transition)
        sapling_mortality = 1 .- sapling_survival
        growth_forms = demography_df.growth_form
        g_jabowas = demography_df.g_jabowa
        max_heights = demography_df.max_hgt
        max_dbhs = demography_df.max_dbh 
        max_ages = demography_df.max_age

        edge_responses = demography_df.edge_response

        #*Based on equations from Botkin et al. (1972)
        b2_jabowas = (2 .* (max_heights .- 1.37)) ./ max_dbhs 
        b3_jabowas = ((max_heights .- 1.37) ./ (max_dbhs).^2) 

        repro_ages = demography_df.repro_age
        repro_heights = demography_df.repro_height
        seed_prod = demography_df.seed_prod
        ldd_dispersal_fracs = demography_df.ldd_dispersal_frac
        ldd_dispersal_dist = demography_df.ldd_dispersal_dist

        regen_heights = demography_df.regen_height

        seedling_inhibition = demography_df.seedling_inhibition

        edge_b0 = zero(Int64)
        edge_b1 = one(Int64)

        seed_list = seeds = Int64[]
        sap_list = Int64[]
        for i in demography_df.growth_form
            seeds = i == 1 ? 10 : 6
            saps = i == 1 ? 2 : 1

            push!(seed_list, seeds)
            push!(sap_list, saps)
        end

        seed_density = sum(seed_list)
        sap_density = sum(sap_list)

        max_shell = max_shade_distance / cell_grain

        n_species = max(n_species, length(seed_prod))
        external_species = demography_df.external_species

        herbivory_amount = demography_df.herbivory

        base_mortality = (4 ./ max_ages)
        supp_tolerance = demography_df.supp_tolerance
        supp_mortality = demography_df.supp_mortality
        gap_maker = demography_df.gap_maker

        shade_tolerance = demography_df.shade_tolerance

        saplings_to_plant = Int64[1,0,1,0,0,1,0,0] #TODO: Make this a parameter or something NB

        #* Is the species able to get soil disease?
        phytothera_target = demography_df.susceptible_soil_disease::Vector{Int64}

        ###--------------------ASSIGN INITIAL PROPERTY VARIABLES----------------------###
        properties = Dict(
            #% PATCH VARIABLES----------------------------#
            :patch_ID => zeros(Int64, prod((dims, dims))),
            :pcor => fill(Tuple{Int64, Int64}[], prod((dims, dims))),
            :seedlings => fill(Int64[], prod((dims, dims))),
            :saplings => fill(Int64[], prod((dims, dims))),
            :edge_weight => zeros(Float64, prod((dims, dims))),
            :previous_species => zeros(Float64, prod((dims, dims))),
            :previous_height => zeros(Float64, prod((dims, dims))),
            :nhb_set => fill(Tuple{Int64, Int64}[], prod((dims, dims))),
            :nhb_set_ids => fill(Int64[], prod((dims, dims))),
            :close_nhbs_count => zeros(Int, prod((dims, dims))), #rename of netlogo models nhbs which is a count of the nearest layer of neighbours
            :nhb_shade_height => zeros(Float64, prod((dims, dims))),
            :nhb_light => zeros(Float64, prod((dims, dims))),
            :disturbed => falses(prod((dims, dims))),
            :expand => falses(prod((dims, dims))),
            :last_change_tick => zeros(Int64, prod((dims, dims))),
            :n_changes => zeros(Int64, prod((dims, dims))),
            :seedling_density => fill(seed_density, prod((dims, dims))), #Could maybe be remvoed and made a reporter using seedlings 
            :sapling_density => fill(sap_density, prod((dims, dims))), #Same as above
            #% GLOBAL VARIABLES---------------------------#
            :tick => zero(1),
            :n_species => n_species::Int64,
            :seedling_survival => seedling_survival::Vector{Float64},
            :sapling_survival => sapling_survival::Vector{Float64},
            :seedling_transition => seedling_transition::Vector{Float64},
            :seedling_mortality => seedling_mortality::Vector{Float64},
            :sapling_mortality => sapling_mortality::Vector{Float64},
            :seedling_inhibition => seedling_inhibition::Vector{Int64},
            :abundances => zeros(Int64, n_species),
            :edge_b0 => edge_b0::Int64,
            :edge_b1 => edge_b1::Int64,
            :growth_forms => growth_forms::Vector{Int64},
            :g_jabowas => g_jabowas::Vector{Float64},
            :max_heights => max_heights::Vector{Int64},
            :max_dbhs => max_dbhs::Vector{Float64},
            :max_ages => max_ages::Vector{Int64},
            :shell_layers => Int64(max_shell),
            :shell_layers_count => collect(range(0, 32, step = 4))::Vector{Int64},
            :base_mortality => base_mortality::Vector{Float64},
            :b2_jabowas => b2_jabowas::Vector{Float64},
            :b3_jabowas => b3_jabowas::Vector{Float64},
            :repro_ages => repro_ages::Vector{Int64},
            :repro_heights => repro_heights::Vector{Float64},
            :seed_prod => seed_prod::Vector{Int64},
            :ldd_dispersal_fracs => ldd_dispersal_fracs::Vector{Float64},
            :ldd_dispersal_dist => ldd_dispersal_dist::Vector{Int64},
            :regen_heights => regen_heights::Vector{Int64},
            :external_species => external_species::Vector{Float64},
            :herbivory_amount => herbivory_amount::Vector{Float64},
            :supp_tolerance => supp_tolerance::Vector{Float64},
            :supp_mortality => supp_mortality::Vector{Float64},
            :gap_maker => gap_maker::Vector{Int64},
            :shade_tolerance => shade_tolerance::Vector{Float64},
            :saplings_to_plant => saplings_to_plant::Vector{Int64},
            :max_density => sap_density::Int64,
            :phytothera_target => phytothera_target::Vector{Int64},
            #% USER INPUTS--------------------------------#
            :cell_grain => cell_grain::Int64,
            :disturbance_freq => disturb_freq::Float64,
            :max_disturb_size => max_disturb_size::Float64,
            :comp_multiplier => comp_multiplier::Float64,
            :edge_effects => edge_effects::Bool,
            :edge_responses => edge_responses::Vector{Float64},
            :external_rain => external_rain::Bool,
            :ext_dispersal_scenario => ext_dispersal_scenario::String,
            :herbivory => herbivory::Bool,
            :saplings_eaten => saplings_eaten::Bool,
            :macro_litter_effect => macro_litter_effect::Float64,
            :ddm => ddm::Bool,
            :restoration_planting => restoration_planting::Bool,
            :planting_frequency => planting_frequency::Int64,
            :phytothera => phytothera::Bool,
            :phyto_global_prob => phyto_global_infection_prob::Float64,
            :phyto_local_prob => phyto_local_infection_prob::Float64,
            :phyto_infectious_radius => phyto_infectious_radius::Int64,
            :phyto_symptoms_dev_prob => phyto_symptoms_dev_prob::Float64,
            :phyto_mortality_prob => phyto_mortality_prob::Float64,
        )

        ###------------------------------CREATE THE MODEL-----------------------------###
        #! Note the type of scheduler can have a large impact on the speed and behaviour
        #! of the final model
        model = ABM(Tree, space; 
            properties,
            rng,
            scheduler = Schedulers.fastest)

        grid = collect(positions(model))
        num_positions = prod((dims, dims))

        ###-----------------UPDATE INITIAL PATCH VALUES TO TRUE VALUES----------------###
        #*For loop runs across all patches and calculates the correct values for each
        #* patch variable
        for p in 1:num_positions
            model.patch_ID[p] = p
            model.pcor[p] = grid[[p]]

            #*Calculate the number of seedlings and saplings for each patch
            model.seedlings[p] = copy(seed_list)
            model.saplings[p] = copy(sap_list)

            #* Calculate the species ID based on the species initial abundance in site_df
            #? Could we use dictionary keys to get name value pairs and make it clearer what we are doing
            #! Column 1 is species column 2 is initial abundance
            specID = wsample(site_df[ : , 1], site_df[ : , 2])

            grow_form = demography_df.growth_form[specID]
            phytothera_infected = false
            phytothera_symptomatic = false


            #% ADD A SINGLE UNIQUE AGENT TO THE PATCH-------------------------#
            #*Use custom function to generate agent dbh, age, and height
            agent_demog = assign_demographic(model,
                                             specID, 
                                             site_df)
            
            add_agent!(grid[p]::Tuple{Int64, Int64}, 
                model, 
                specID::Int64, 
                p::Int64, #patch_here_ID,
                grow_form::Int64, 
                agent_demog[1]::Float64, #height
                agent_demog[2]::Float64, #dbh
                agent_demog[3]::Float64, #age
                Float64[],
                phytothera_infected::Bool,
                zero(Int64), #phytothera_infected_age
                phytothera_symptomatic::Bool,
                )

            #% UPDATE PATCH LEVEL PROPERTIES----------------------------------#
            #* calculate edge distance
            e_dist = minimum(grid[p] .- minimum(positions(model)))
            #! Note I have renamed edge-b2 from Netlogo to edge_strength
            weight = edge_b1 .* exp(-edge_strength .* e_dist) .+ edge_b0
            model.edge_weight[p] = weight

            model.nhb_set[p] = collect(nearby_positions(grid[p], 
                                                        model::ABM{<:GridSpaceSingle}, 
                                                        Int64(max_shell)))

            model.close_nhbs_count[p] = length(collect(nearby_positions(grid[p], 
                                                                        model::ABM{<:GridSpaceSingle}, 
                                                                        1)))
        end

        ###----------------------------UPDATE SECOND STAGE VARIABLES----------------------------###
        #! Note we have a second loop to assign some features as by default things in Julia do not
        #! seem to run sequentially but rather all together, hence the ability to assign a function
        #! after calling it but this means trying to get nhb_set in the same loop they are assigned 
        #! does not seem to function.
        pcors = model.pcor
        nhb_sets = model.nhb_set
        crit_heights = range(0, 32, step = 4)

        Threads.@threads for i in 1:num_positions
            model.nhb_shade_height[i] = set_get_functions.get_nhb_shade_height(i, 
                                                                               model,
                                                                               grid,
                                                                               crit_heights,
                                                                               Int64(max_shell))

            n_ids = Int64[]
            for n in nhb_sets[i]
                n_id = findfirst(isequal([n]), pcors)
                push!(n_ids, n_id)
            end
            model.nhb_set_ids[i] = n_ids
        end
        
        return model
    end


    #//-------------------------------------------------------------------------------------------#
    #% CALCULATE AGENT DEMOGRAPHIC VALUES FUNCTION
    """
    # Report agent demographic parameters
    Function calculates an initial diameter at breast height, height and age for an agent based
    on the agents species.
    ## Arguments:
    - `model`: The AgentBasedModel object defining the current model. This object is usually 
    created using `Agents.ABM()`.
    - `species::Integer`: Selected species ID.
    - `site_df::DataFrame`: Data defining the initial site properties to be modelled
    """
    function assign_demographic(
        model,
        species::Integer,
        site_df::DataFrame = site_df
    )
        ###-----------------DEFINE SPECIES CHARACTERISTICS------------------###
        growth_form = model.growth_forms[species]

        max_height_frac = site_df.max_init_hgt[species]
        max_height_frac = max_height_frac < 0 || max_height_frac > 1 ? 0.95 : max_height_frac
        max_height = model.max_heights[species]

        max_dbh = model.max_dbhs[species]
        start_dbh = site_df.start_dbh[species]
        start_dbh_sd = site_df.start_dbh_sd[species]

        b2_jabowa = model.b2_jabowas[species]
        b3_jabowa = model.b3_jabowas[species]
        g_jabowa = model.g_jabowas[species]

        #% DEFINE BEHAVIOUR FOR TREES (GROWTH FORM 1)---------------#
        if growth_form == 1
            ## Define initial DBH
            dbh = min(rand(distribution_functions.generate_LogNormal(start_dbh,
                                                                     start_dbh_sd), 1)[1], 
                (max_height_frac * max_dbh))
            dbh = max(0.01, dbh)

            ## Define initial height
            height = 1.37 + (b2_jabowa * dbh) - (b3_jabowa * dbh * dbh)

            ## Define initial age
            age = demog_metrics.age_by_dbh(
                height, 
                dbh,
                max_dbh,
                Float64(max_height),
                g_jabowa,
                b2_jabowa,
                b3_jabowa
            )


        #% DEFINE BEHAVIOUR FOR TREE FERNS (GROWTH FORM 2)----------#
        elseif growth_form == 2
            ## Define initial height
            height = min(rand(distribution_functions.generate_LogNormal(start_dbh,
                                                                        start_dbh_sd), 1)[1], 
                    (max_height_frac * max_height))
            
            #* Ensure height is non-negative    
            height = height < 0 ? 1.5 : height

            ## Define initial DBH (Note this is a staitc value)
            dbh = 0.1

            ## Define initial age
            age = demog_metrics.age_by_height(
                height
            )

        # If growth form is not of known type send an error message
        else
            error("The growth form $growth_form is undefined, please check species demography data")
        end

        #! Note order of output is important as it is used implicitly throughout model
        return(height, dbh, age)
    end
end