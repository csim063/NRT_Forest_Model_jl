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
        growth_form::Int
        height::Float64
        dbh::Float64
        age::Float64
    end

    #% Define world
    function forest_model(;
        forest_area = 16,
        cell_grain = 4,
        edge_strength = 0.0,
        max_shade_distance = 32,
        site_df = site_df,
        demography_df = demography_df,
        seed = 999,)
        
        #? Could we define space as many dimensional and just add the properties like height
        #? and species present to this
        ## Define globals
        dims = trunc(Int, (sqrt(forest_area * 1e4)) / cell_grain)

        space = GridSpaceSingle((dims, dims); periodic = false, metric = :manhattan);
        rng = MersenneTwister(seed)

        seedling_survival = demography_df.seedling_survival
        sapling_survival = demography_df.sapling_survival
        seedling_transition = demography_df.seedling_transition
        seedling_mortality = 1 .- (seedling_survival .+ seedling_transition)
        sapling_mortality = 1 .- sapling_survival

        edge_b0 = 0
        edge_b1 = 1

        seed_list = Int64[]
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

        ## Define patch properties
        properties = (
            seedlings = fill(seed_list, prod((dims, dims))),
            saplings = fill(sap_list, prod((dims, dims))),
            edge_weight = zeros(Float64, prod((dims, dims))),
            previous_species = fill(Int64[], prod((dims, dims))),
            previous_height = fill(Int64[], prod((dims, dims))),
            previous_growth = fill(Float64[], prod((dims, dims))),
            nhb_set = fill(fill(Tuple{Int64, Int64}[], Int(max_shell)), prod((dims, dims))),
            close_nhbs_count = zeros(Int, prod((dims, dims))), #rename of netlogo models nhbs which is a count of the nearest layer of neighbours
            nhb_shade_height = zeros(Float64, prod((dims, dims))),
            nhb_light = zeros(Float64, prod((dims, dims))),
            disturbed = zeros(Float64, prod((dims, dims))),
            expand = falses(prod((dims, dims))),
            last_change_tick = fill(Int64[], prod((dims, dims))),
            n_changes = zeros(Float64, prod((dims, dims))),
            seedling_density = fill(seed_density, prod((dims, dims))), #Could maybe be remvoed and made a reporter using seedlings 
            sapling_density = fill(sap_density, prod((dims, dims))) #Same as above
        )
        model = ABM(Tree, space; 
            properties,
            rng,
            scheduler = Schedulers.Randomly())

        ## Populate the world with adult tree agents
        grid = collect(positions(model))
        num_positions = prod((dims, dims))

        #Make for loop that samples a proportion of space and allocates each species
        for p in 1:num_positions
            #? Could we use dictionary keys to get name value pairs and make it clearer what we are doing
            # Column 1 is species column 2 is initial abundance
            specID = wsample(site_df[ : , 1], site_df[ : , 2])

            grow_form = demography_df.growth_form[specID]

            ## Get height dbh and age
            agent_demog = assign_demographic(specID, 
                                             site_df, 
                                             demography_df)

            adult_tree = Tree(
                p,
                grid[p],
                specID,
                grow_form,
                agent_demog[1],
                agent_demog[2],
                agent_demog[3],
            )
            add_agent_single!(adult_tree, model)

            ## Update patch level properties
            e_dist = minimum(grid[p] .- minimum(positions(model)))
            weight = edge_b1 * exp(-edge_strength * e_dist) + edge_b0
            model.edge_weight[p] = weight

            #TODO fix this to be vectorised
            d_nhbs = fill(Tuple{Int64, Int64}[], Int(max_shell))
            for d in range(1, Int(max_shell))
                nhbs_ids = Tuple{Int64, Int64}[]
                for idx in nearby_positions(grid[p], model::ABM{<:GridSpaceSingle}, d)
                    push!(nhbs_ids, idx)
                end
                d_nhbs[d] = d ≤ 1 ? nhbs_ids : setdiff(nhbs_ids, d_nhbs[d - 1])
            end
            model.nhb_set[p] = d_nhbs
            model.close_nhbs_count[p] = length(d_nhbs[1])
        end

        #! Note we have a second loop to assign some features as by default things in Julia do not
        #! seem to run sequentially but rather all together, hence the ability to assign a function
        #! after calling it but this means trying to get nhb_set in the same loop they are assigned 
        #! does not seem to function.
        for i in 1:num_positions
            model.nhb_shade_height[i] = set_get_functions.get_nhb_shade_height(i, model)
            #model.nhb_light[i] = set_get_functions.get_light_env(i, model, demography_df) #* Not used in initialisation
        end
        
        return model
    end

    #// Define patches

    function assign_demographic(
        species::Integer,
        site_df = site_df,
        demography_df = demography_df
    )
        ## Define species charatecteristics
        growth_form = demography_df.growth_form[species]

        max_height_frac = site_df.max_init_hgt[species]
        max_height_frac = max_height_frac < 0 || max_height_frac > 1 ? 0.95 : max_height_frac
        max_height = demography_df.max_hgt[species]

        max_dbh = demography_df.max_dbh[species]
        start_dbh = site_df.start_dbh[species]
        start_dbh_sd = site_df.start_dbh_sd[species]

        b2_jabowa = (2 * (max_height - 1.37)) / max_dbh #*Based on equation from Botkin 2001
        b3_jabowa = ((max_height - 1.37) / (max_dbh)^2) #*Based on equation from Botkin 2001
        g_jabowa = demography_df.g_jabowa[species]

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
