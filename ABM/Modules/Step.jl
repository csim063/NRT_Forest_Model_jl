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
        nhb_shade_height = model.nhb_shade_height,
        comp_multiplier = model.comp_multiplier,
        edge_effects = model.edge_effects,
        edge_weights = model.edge_weight,
        edge_responses = model.edge_responses,
        g_jabowas = model.g_jabowas,
        b2_jabowas = model.b2_jabowas,
        b3_jabowas = model.b3_jabowas,
        max_dbhs = model.max_dbhs,
        max_heights = model.max_heights
    )
    
        #% DEFINE VARIABLES USED ACROSS PROCEDURES------------------#
        #* These variables are defined in the agent stepping
        #* function itself as they are unique to each agent.
        spec_num = agent.species_ID
        cell = agent.patch_here_ID
        age = agent.age
        shade_height = nhb_shade_height[spec_num]
        edge_weight = edge_weights[spec_num]
        edge_response = edge_responses[spec_num]
        g_jabowa = g_jabowas[spec_num]
        b2_jabowa = b2_jabowas[spec_num]
        b3_jabowa = b3_jabowas[spec_num]
        max_dbh = max_dbhs[spec_num]
        max_height = max_heights[spec_num]

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
            if sp > zero(1) #*zero(1) gives a 0 value which is more type stable than 0
                nhbs_ids = model.nhb_set_ids[agent.patch_here_ID]

                demog_funcs.nhb_dispersal(model,
                                        sp,
                                        ldd_disp_frac,
                                        r_hgt,
                                        agent.dbh,
                                        cell_grain,
                                        model.shell_layers,
                                        a_position, 
                                        nhbs_ids,
                                        shad_hs,
                                        seedlings,
                                        spec_num
                                        )

            end

            ldd_disp_dist = model.ldd_dispersal_dist[agent.species_ID]
            demog_funcs.ldd_within(model,
                                   sp,
                                   ldd_disp_frac,
                                   ldd_disp_dist,
                                   cell_grain,
                                   a_position,
                                   p_cors,
                                   shad_hs,
                                   r_hgt,
                                   seedlings,
                                   spec_num)
        end

        #% MORTAILTY------------------------------------------------#
        ## Herbivory
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
                                model.gap_maker,
                                model.expand,
                                model.previous_species,
                                model.previous_height,
                                agent.height,
                                agent.id,
                                age,
                                agent.previous_growth,
                                model.supp_tolerance,
                                model.supp_mortality
                                )
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
                                                  saplings
                                                  )
        end

        #% BEYOND PATCH DISPERSAL-----------------------------------#
        if external_rain == true
            demog_funcs.external_ldd(ext_dispersal_scenario,
                                     grid,
                                     n_species,
                                     abundances,
                                     external_species,
                                     seedlings)
        end

        #% EXPAND FOREST GAPS---------------------------------------#
        #* Only current forest gap patches may check if the gap may 
        #* expand
        for ep in model.patch_ID[model.expand .== true]
            demog_funcs.expand_gap(ep,
                                   model,
                                   grid)    

            model.expand[ep] = false
        end

        for i in 1:length(grid)
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
        for _ in eachindex(empty_patches)
            patch = random_empty(model)
            cell_ID = findfirst(isequal([patch]), pcor)

            demog_funcs.capture_gap([cell_ID], 
                                    seedlings,
                                    saplings,
                                    nhb_light,
                                    shade_tolerance,
                                    growth_forms,
                                    b2_jabowas,
                                    b3_jabowas,
                                    model.last_change_tick,
                                    tick,
                                    model.n_changes,
                                    new_agents_list
                                    )

        end

        #* capture_gap() only adds the properties of a new agent to assign
        #* to a list this next loop actually assigns those agents to patches.
        empty_patches = Random.shuffle!(empty_patches)
        for e in eachindex(new_agents_list)
            pos = empty_patches[e]
            agent = new_agents_list[e]

            add_agent!(
                    pos,
                    model,
                    agent[1], 
                    agent[2],
                    agent[3], 
                    agent[4], 
                    agent[5], 
                    agent[6], 
                    Float64[]
                    )
        end

        #% RESTORATION----------------------------------------------#
        if restoration_planting == true && mod(tick, planting_frequency) == 0
            for i in 1:length(grid)
                saplings[i] .+= saplings_to_plant
            end
        end

        #% SEEDLING AND SAPLING MORTAILTY AND GROWTH----------------#
        for i in 1:length(grid)
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

        model.tick += 1
    end
end