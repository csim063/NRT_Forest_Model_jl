"""
Module contains the agent and model stepping functions which are the key functions which actually
run the model. Can be thought of as the "go" procedure from NetLogo.
"""

module go
    using Agents
    using Random
    using StatsBase
    using DataFrames
    using Distributions

    #* CUSTOM MODULES
    include("Disturbance_functions.jl")
    include("Demographic_functions.jl")
    include("Helper_functions.jl")
    include("Disease_functions.jl")
    
    #//-------------------------------------------------------------------------------------------#
    #% AGENT STEPPING FUNCTION
    """
    This function applies procedures to each agent/tree in the model one by one. Note that as each
    agent stepping function is applied to each agent it is run multiple times per tick if there are
    multiple agents.
    """
    function agent_step!(
        agent,
        model,
        nhb_shade_height = model.nhb_shade_height::Vector{Float64},
        comp_multiplier = model.comp_multiplier::Float64,
        edge_effects = model.edge_effects::Bool,
        edge_weights = model.edge_weight::Vector{Float64},
        edge_responses = model.edge_responses::Vector{Float64},
        g_jabowas = model.g_jabowas::Vector{Float64},
        b2_jabowas = model.b2_jabowas::Vector{Float64},
        b3_jabowas = model.b3_jabowas::Vector{Float64},
        max_dbhs = model.max_dbhs::Vector{Float64},
        max_heights = model.max_heights::Vector{Int64},
        expand = model.expand::BitVector,
        phytothera = model.phytothera::Bool,
        phytothera_target = model.phytothera_target::Vector{Int64},
        phyto_global_prob = model.phyto_global_prob::Float64,
        phyto_local_prob = model.phyto_local_prob::Float64,
        phyto_infectious_radius = model.phyto_infectious_radius::Int64,
        phyto_symptoms_dev_prob = model.phyto_symptoms_dev_prob::Float64,
        phyto_mortality_prob = model.phyto_mortality_prob::Float64,
        phyto_transmission_age = model.phyto_transmission_age::Int64,
        phyto_min_symptomatic_age = model.phyto_min_symptomatic_age::Int64,
        rust = model.rust::Bool,
        rust_target = model.rust_target::Vector{Int64},
        rust_infection_prob = model.rust_global_infection_prob::Float64,
        rust_symptoms_dev_prob = model.rust_symptoms_dev_prob::Float64,
        rust_mortality_prob = model.rust_mortality_prob::Float64,
        rust_min_symptomatic_age = model.rust_min_symptomatic_age::Int64,
        grass::Bool = model.grass,
        grass_invasion_prob::Float64 = model.grass_invasion_prob,
        grass_colonisation_prob::Float64 = model.grass_colonisation_prob,
    )
    
        #% DEFINE VARIABLES USED ACROSS PROCEDURES------------------#
        #* These variables are defined in the agent stepping
        #* function itself as they are unique to each agent.
        spec_num::Int64 = agent.species_ID
        cell::Int64 = agent.patch_here_ID
        age::Float64 = agent.age
        tree_ID::Int64 = agent.id
        agent_pos::Tuple{Int64, Int64} = agent.pos

        shade_height::Float64 = nhb_shade_height[spec_num]
        edge_weight::Float64 = edge_weights[spec_num]
        edge_response::Float64 = edge_responses[spec_num]
        g_jabowa::Float64 = g_jabowas[spec_num]
        b2_jabowa::Float64 = b2_jabowas[spec_num]
        b3_jabowa::Float64 = b3_jabowas[spec_num]
        max_dbh::Float64 = max_dbhs[spec_num]
        max_height::Int64 = max_heights[spec_num]
        gap_maker::Int64 = model.gap_maker[spec_num]
        pest_herbivory::Float64 = model.adult_pest_herbivory[spec_num]

        #% GROW-----------------------------------------------------#
        #*Have each tree grow, i.e. increase their age, height and 
        #* diameter at breast height.
        demog_funcs.grow(agent, 
                        shade_height,
                        agent.height,
                        comp_multiplier,
                        edge_effects,
                        edge_weight,
                        edge_response,
                        agent.growth_form,
                        agent.dbh,
                        g_jabowa,
                        b2_jabowa,
                        b3_jabowa,
                        max_dbh,
                        max_height)

        #% DISPERSAL------------------------------------------------#
        #* This is only dispersal of seeds local to the habitat
        #* patch not dispersal rain from beyond the patch
        if age ≥ model.repro_ages[spec_num] && agent.height ≥ model.repro_heights[spec_num]
            sp = Int64(model.seed_prod[agent.species_ID])
            ldd_disp_frac = model.ldd_dispersal_fracs[agent.species_ID]
            cell_grain = model.cell_grain
            a_position = agent.pos
            p_cors = model.pcor
            shad_hs = nhb_shade_height
            r_hgt = model.regen_heights[agent.species_ID]
            seedlings = model.seedlings
            if (grass == true && rand() < grass_colonisation_prob) || grass == false
                if sp > zero(1) #*zero(1) gives a 0 value which is more type stable than 0
                    nhbs_ids = model.nhb_set_ids[agent.patch_here_ID]

                    demog_funcs.nhb_dispersal(model,
                                            sp::Int64,
                                            ldd_disp_frac::Float64,
                                            r_hgt::Int64,
                                            agent.dbh::Float64,
                                            cell_grain::Int64,
                                            model.shell_layers::Int64,
                                            a_position::Tuple{Int64, Int64},
                                            nhbs_ids::Vector{Int64},
                                            shad_hs::Vector{Float64},
                                            seedlings::Vector{Vector{Int64}},
                                            spec_num::Int64
                                            )

                end
            end

            if (grass == true && rand() < grass_colonisation_prob) || grass == false
                ldd_disp_dist = model.ldd_dispersal_dist[agent.species_ID]
                demog_funcs.ldd_within(model,
                                    sp::Int64,
                                    ldd_disp_frac::Float64,
                                    ldd_disp_dist::Int64,
                                    cell_grain::Int64,
                                    a_position::Tuple{Int64,Int64},
                                    p_cors::Vector{Vector{Tuple{Int64, Int64}}},
                                    shad_hs::Vector{Float64},
                                    r_hgt::Int64,
                                    seedlings::Vector{Vector{Int64}},
                                    spec_num::Int64)
            end
        end

        #% MORTAILTY------------------------------------------------#
        ## Herbivory for seedlings and saplings
        if model.herbivory == true
            demog_funcs.herbivore_effect(agent, model)
        end

        ## Seedling and sapling background mortality and growth
        demog_funcs.thin_regenbank(agent, model)

        ## Mortality due to macro-litterfall
        #* This is only done under tree-ferns (growth form 2)
        if agent.growth_form == 2
            demog_funcs.macro_litter_fall(agent, model)
        end

        ## Background and competition mortality
        #* Due to split between agent and model stepping functions
        #* a check to ensure that agents have not just been 
        #* created (i.e. age ≥ 1) is needed.
        if agent.age ≥ 1.0
            demog_funcs.death(agent, 
                                model,
                                cell,
                                model.ddm,
                                spec_num,
                                model.base_mortality,
                                gap_maker,
                                expand,
                                model.previous_height,
                                agent.height,
                                tree_ID,
                                age,
                                agent.previous_growth,
                                model.supp_tolerance,
                                model.supp_mortality,
                                grass,
                                grass_invasion_prob,
                                model.herbivory,
                                pest_herbivory,
                                model.herbivore_variability,
                                model.weather_adjustment
                                )
        end

        #% DISEASE--------------------------------------------------#
        ## Phytothera (e.g. Kauri dieback)
        if phytothera == true && id_in_position(agent_pos, model::ABM{<:GridSpaceSingle}) != 0
            phytothera_infected::Bool = agent.phytothera_infected

            #* Spread the infection to non-infected target trees
            if phytothera_target[spec_num] == 1 && phytothera_infected == false
                disease_functions.phytothera_spread(agent,
                                                    model,
                                                    phyto_global_prob::Float64,
                                                    phyto_local_prob::Float64,
                                                    phyto_infectious_radius::Int64,
                                                    phyto_transmission_age::Int64,
                                                    agent_pos::Tuple{Int64, Int64},
                                                    edge_effects::Bool,
                                                    )
            end

            #* Apply the disease effects to infected trees
            if phytothera_infected == true
                disease_functions.phytothera_impact(agent,
                                                    model,
                                                    agent.phytothera_infected::Bool,
                                                    phyto_min_symptomatic_age::Int64,
                                                    phyto_symptoms_dev_prob::Float64,
                                                    phyto_mortality_prob,
                                                    gap_maker::Int64,
                                                    spec_num::Int64,
                                                    cell::Int64,
                                                    expand::BitVector,
                                                    model.previous_height::Vector{Float64},
                                                    agent.height::Float64,
                                                    tree_ID::Int64,
                                                    grass::Bool,
                                                    grass_invasion_prob::Float64,)
            end
        end 

        ## Rust (e.g. Myrtle rust)
        if rust == true && id_in_position(agent_pos, model::ABM{<:GridSpaceSingle}) != 0
            rust_infected::Bool = agent.rust_infected

            #* Spread the infection to non-infected target trees
            if rust_target[spec_num] == 1 && rust_infected == false
                disease_functions.rust_spread(agent,
                                              rust_infection_prob)
            end

            #* Apply the disease effects to infected trees
            if rust_infected == true
                disease_functions.rust_impact(agent,
                                              rust_infected,
                                              rust_min_symptomatic_age,
                                              rust_symptoms_dev_prob,
                                              rust_mortality_prob,
                                              gap_maker,
                                              spec_num,
                                              cell,
                                              expand,
                                              model.previous_height,
                                              agent.height,
                                              tree_ID,
                                              grass,
                                              grass_invasion_prob,)
            end
        end
    end

    #//-------------------------------------------------------------------------------------------#
    #% MODEL STEPPING FUNCTION
    #* RUN ONCE PER MODEL PER CALL (I.E. ONCE PER TICK)
    """
    This function applies procedures that are not applicable to be run by a single agent, i.e.
    procedures that affect patches and globals. This function is run strictly once per model tick.
    By default the model stepping function is run after the agent stepping function.
    """
    function model_step!(
                        model,
                        tick::Int64 = model.tick,
                        saplings::Vector{Vector{Int64}} = model.saplings,
                        seedlings::Vector{Vector{Int64}} = model.seedlings,
                        sapling_density::Vector{Int64} = model.sapling_density,
                        seedling_density::Vector{Int64} = model.seedling_density,
                        n_species::Int64 = model.n_species,
                        disturbance_freq::Float64 = model.disturbance_freq,
                        disturbed::BitVector = model.disturbed,
                        max_disturb_size::Float64 = model.max_disturb_size,
                        nhb_set_ids::Vector{Vector{Int64}} = model.nhb_set_ids,
                        close_nhbs_count::Vector{Int64} = model.close_nhbs_count,
                        patch_ID::Vector{Int64} = model.patch_ID,
                        shell_layer_count::Vector{Int64} = model.shell_layers_count,
                        shell_layers::Int64 = model.shell_layers,
                        external_rain::Bool = model.external_rain,
                        ext_dispersal_scenario::String = model.ext_dispersal_scenario,
                        abundances::Vector{Int64} = model.abundances,
                        external_species::Vector{Float64} = model.external_species,
                        max_heights::Vector{Int64} = model.max_heights,
                        pcor::Vector{Vector{Tuple{Int64, Int64}}} = model.pcor,
                        shade_tolerance::Vector{Float64} = model.shade_tolerance,
                        growth_forms::Vector{Int64} = model.growth_forms,
                        b2_jabowas::Vector{Float64} = model.b2_jabowas,
                        b3_jabowas::Vector{Float64} = model.b3_jabowas,
                        saplings_to_plant::Vector{Int64} = model.saplings_to_plant,
                        restoration_planting::Bool = model.restoration_planting,
                        planting_frequency::Int64 = model.planting_frequency,
                        seedling_mortality::Vector{Float64} = model.seedling_mortality,
                        sapling_mortality::Vector{Float64} = model.sapling_mortality,
                        seedling_transition::Vector{Float64} = model.seedling_transition,
                        grass::Bool = model.grass,
                        grass_invasion_prob::Float64 = model.grass_invasion_prob,
                        grass_colonisation_prob::Float64 = model.grass_colonisation_prob,
                        weather::Bool = model.weather,
                        weather_variability::Float64 = model.weather_variability,
                        weather_adjustment::Float64 = model.weather_adjustment,
                        )
        #% DEFINE VARIABLES USED ACROSS PROCEDURES------------------#
        grid = collect(positions(model))

        #% DISTURBANCE EVENTS---------------------------------------#
        if model.disturbance_freq > 0
            disturbance_functions.lsp_disturbance(model,
                                                  grid,
                                                  disturbance_freq,
                                                  disturbed,
                                                  max_disturb_size,
                                                  nhb_set_ids,
                                                  close_nhbs_count,
                                                  patch_ID,
                                                  n_species,
                                                  seedlings,
                                                  saplings,
                                                  grass,
                                                  grass_invasion_prob
                                                  )
        end

        #% BEYOND PATCH DISPERSAL-----------------------------------#
        if external_rain == true
            if (grass == true && rand() < grass_colonisation_prob) || grass == false
                demog_funcs.external_ldd(ext_dispersal_scenario,
                                        grid,
                                        n_species,
                                        abundances,
                                        external_species,
                                        seedlings,)
            end
        end

        #% EXPAND FOREST GAPS---------------------------------------#
        #* Only current forest gap patches may check if the gap may 
        #* expand
        for ep in model.patch_ID[model.expand .== true]
            demog_funcs.expand_gap(ep,
                                   model,
                                   grid,
                                   grass,
                                   grass_invasion_prob)    

            model.expand[ep] = false
        end

        Threads.@threads for i in eachindex(grid)
            #% NEIGHBOURHOOD FUNCTIONS------------------------------#
            ##TODO DOES THIS NEED TO BE DONE FOR EVERY PATCH?
            #* Use the neighbours sets for all cells to determine
            #* there mean shade height and light environment.
            model.nhb_shade_height[i] = set_get_functions.get_nhb_shade_height(i, 
                                                                               model,
                                                                               grid,
                                                                               shell_layer_count,
                                                                               shell_layers)

            model.nhb_light[i] = set_get_functions.get_light_env(model.nhb_shade_height[i], 
                                                                max_heights)
        end

        #% GROW NEW TREES IN GAPS-----------------------------------#
        #TODO: This is a bit of a mess, needs to be tidied up and made into a function
        #* Define pre-loop variables
        empty_patches = Random.shuffle!(collect(empty_positions(model))) 
        nhb_light = model.nhb_light
        new_agents_list = Any[]

        #* Loop through all empty patches and grow trees in gaps
        for _ in eachindex(empty_patches::Vector{Tuple{Int64, Int64}})
            patch = random_empty(model)
            cell_ID = findfirst(isequal([patch::Tuple{Int64, Int64}]), 
                                        pcor::Vector{Vector{Tuple{Int64, Int64}}}
                                        )

            demog_funcs.capture_gap([cell_ID::Int64], 
                                    seedlings::Vector{Vector{Int64}},
                                    saplings::Vector{Vector{Int64}},
                                    nhb_light::Vector{Float64},
                                    shade_tolerance::Vector{Float64},
                                    growth_forms::Vector{Int64},
                                    b2_jabowas::Vector{Float64},
                                    b3_jabowas::Vector{Float64},
                                    model.last_change_tick::Vector{Int64},
                                    tick::Int64,
                                    model.n_changes::Vector{Int64},
                                    new_agents_list::Vector{Any},
                                    )

        end

        #* capture_gap() only adds the properties of a new agent to assign
        #* to a list this next loop actually assigns those agents to patches.
        empty_patches = Random.shuffle!(empty_patches)
        for e in eachindex(new_agents_list)
            pos = empty_patches[e]
            agent = new_agents_list[e]

            add_agent!(
                    pos::Tuple{Int64, Int64},
                    model,
                    trunc(Int, agent[1])::Int64, 
                    trunc(Int, agent[2])::Int64, 
                    trunc(Int, agent[3])::Int64, 
                    agent[4]::Float64, 
                    agent[5]::Float64, 
                    agent[6]::Float64, 
                    Float64[],
                    false, ## Phytothera
                    zero(Int64),
                    false,
                    false, ## Rust 
                    zero(Int64),
                    false,
                    )
        end

        #% RESTORATION----------------------------------------------#
        if restoration_planting == true && mod(tick, planting_frequency) == 0
            for i in eachindex(grid)
                saplings[i] .+= saplings_to_plant
            end
        end

        #% SEEDLING AND SAPLING MORTAILTY AND GROWTH----------------#
        if weather
            seedling_mortality .+= weather_adjustment
            sapling_mortality .+= weather_adjustment
        end
        Threads.@threads for i in eachindex(grid)
            demog_funcs.regenerate_patch_bank(i,
                                            seedlings, 
                                            saplings,
                                            seedling_mortality,
                                            sapling_mortality,
                                            seedling_transition)
                                            
            seedling_density[i] = sum(seedlings[i])
            sapling_density[i] = sum(saplings[i])
        end

        #% REPORTERS------------------------------------------------#
        model.max_density = maximum(sapling_density)
        set_get_functions.update_abundances(model, n_species)

        #% UPDATE END OF TICK VARIABLES-----------------------------#
        #* Redraw weathers impact on mortality if it is included in the model
        if weather
            model.weather_adjustment = rand(Normal(0, weather_variability))
        end
        model.tick += 1
    end
end