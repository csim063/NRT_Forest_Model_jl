"""
Module contains all of the functions required to setup the model at the start
of a run
"""
module Setup
    using Agents
    using Random
    using StatsBase
    using DataFrames

    #% Helper functions
    include("Helper_functions.jl")
    include("Demographic_assignments.jl")

    #% Define agents
    #? Maybe a dict of species key value pairs could be useful
    @agent Tree GridAgent{2} begin 
        species_ID::Int
        patch_here_ID::Int
        growth_form::Int
        height::Float64
        dbh::Float64
        age::Float64
        previous_growth::Array
    end

    #% Define world
    function forest_model(;
        forest_area = 16,
        cell_grain = 4,
        n_species = 8,
        edge_strength = 0.0,
        max_shade_distance = 32,
        site_df = site_df,
        demography_df = demography_df,
        seed = 999,
        disturb_freq = 0.100,
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
        planting_frequency = 10)
        
        #? Could we define space as many dimensional and just add the properties like height
        #? and species present to this
        ## Define globals
        dims = trunc(Int, (sqrt(forest_area * 1e4)) / cell_grain)

        space = GridSpaceSingle((dims, dims); periodic = false, metric = :chebyshev);
        rng = MersenneTwister(seed)

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

        b2_jabowas = (2 * (max_heights .- 1.37)) ./ max_dbhs #*Based on equation from Botkin 2001
        b3_jabowas = ((max_heights .- 1.37) ./ (max_dbhs).^2) #*Based on equation from Botkin 2001

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

        saplings_to_plant = Int64[1,0,1,0,0,1,0,0]

        properties = Dict(
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
            :expand => falses(prod((dims, dims))),
            #%Globals
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
            :new_agents_list => Any[],
            #%User inputs
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
            :planting_frequency => planting_frequency::Int64
        )
        #TODO Give some thought as to the most efficient scheduler to use
        model = ABM(Tree, space; 
            properties,
            rng,
            scheduler = Schedulers.fastest)

        ## Populate the world with adult tree agents
        grid = collect(positions(model))
        num_positions = prod((dims, dims))

        #Make for loop that samples a proportion of space and allocates each species
        for p in 1:num_positions
            model.patch_ID[p] = p
            model.pcor[p] = grid[[p]]

            model.seedlings[p] = copy(seed_list)
            model.saplings[p] = copy(sap_list)

            #? Could we use dictionary keys to get name value pairs and make it clearer what we are doing
            # Column 1 is species column 2 is initial abundance
            specID = wsample(site_df[ : , 1], site_df[ : , 2])
            patch_here_ID = model.patch_ID[p]

            grow_form = demography_df.growth_form[specID]

            ## Get height dbh and age
            agent_demog = assign_demographic(model,
                                             specID, 
                                             site_df)

            add_agent!(grid[p], model, 
                specID, 
                patch_here_ID,
                grow_form, 
                agent_demog[1], #height
                agent_demog[2], #dbh
                agent_demog[3], #age
                Float64[]
                )

            ## Update patch level properties
            e_dist = minimum(grid[p] .- minimum(positions(model)))
            #! Note I have renamed edge-b2 as edge_strength
            weight = edge_b1 * exp(-edge_strength * e_dist) + edge_b0
            model.edge_weight[p] = weight

            model.nhb_set[p] = collect(nearby_positions(grid[p], 
                                                        model::ABM{<:GridSpaceSingle}, 
                                                        Int64(max_shell)))

            model.close_nhbs_count[p] = length(collect(nearby_positions(grid[p], 
                                                                        model::ABM{<:GridSpaceSingle}, 
                                                                        1)))
        end

        #! Note we have a second loop to assign some features as by default things in Julia do not
        #! seem to run sequentially but rather all together, hence the ability to assign a function
        #! after calling it but this means trying to get nhb_set in the same loop they are assigned 
        #! does not seem to function.
        for i in 1:num_positions
            model.nhb_shade_height[i] = set_get_functions.get_nhb_shade_height(i, 
                                                                               model,
                                                                               collect(positions(model)),
                                                                               range(0, 32, step = 4),
                                                                               Int64(max_shell))

            n_ids = Int64[]
            for n in model.nhb_set[i]
                n_id = findfirst(isequal([n]), model.pcor)
                push!(n_ids, n_id)
            end
            model.nhb_set_ids[i] = n_ids
        end
        
        return model
    end

    #// Define patches

    function assign_demographic(
        model,
        species::Integer,
        site_df = site_df
    )
        ## Define species charatecteristics
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

        ## Define behaviour for trees (growth form 1)
        if growth_form == 1
            # Define initial DBH
            dbh = min(rand(distribution_functions.generate_LogNormal(start_dbh,
                                                                     start_dbh_sd), 1)[1], 
                (max_height_frac * max_dbh))
            dbh = max(0.01, dbh)

            # Define initial height
            height = 1.37 + (b2_jabowa * dbh) - (b3_jabowa * dbh * dbh)

            # Define initial age
            age = demog_metrics.age_by_dbh(
                height, 
                dbh,
                max_dbh,
                Float64(max_height),
                g_jabowa,
                b2_jabowa,
                b3_jabowa
            )


        #* Define behaviour for tree ferns (growth form 2)
        elseif growth_form == 2
            #? Is this actually the correct calculation
            height = min(rand(distribution_functions.generate_LogNormal(start_dbh,
                                                                        start_dbh_sd), 1)[1], 
                    (max_height_frac * max_height))
            height = height < 0 ? 1.5 : height

            ## TODO This hard coding seems odd
            dbh = 0.1

            age = demog_metrics.age_by_height(
                height
            )

        #* If growth form is not of known type send an error message
        else
            error("The growth form $growth_form is undefined, please check species demography data")
        end

        return(height, dbh, age)
    end
end