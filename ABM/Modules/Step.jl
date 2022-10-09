"""
Module contains the step (aka go function) for the model
"""

module go
    using Agents
    using Random
    using StatsBase
    using DataFrames
    using Distributions

    include("Disturbance_functions.jl")
    include("Demographic_functions.jl")
    include("Helper_functions.jl")
    
    
    #TODO ORDERS MAY BE ODD WITH AGENT AND MODEL STEPS SPLITTING GO FUNCTION
    #//-------------------------------------------------------------------------------------------#
    #%This is the step function for the individual trees (no globals changed) 
    #* RUN ONCE PER AGENT PER CALL (I.E. Multiple times per tick if multiple agents)
    """
    TODO some proper function documentation
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
            demog_funcs.ldd_within(agent, 
                                   model,
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

        if model.herbivory == true
            demog_funcs.herbivore_effect(agent, model)
        end

        demog_funcs.thin_regenbank(agent, model)

        #* Run for tree ferns only
        if agent.growth_form == 2
            demog_funcs.macro_litter_fall(agent, model)
        end

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


    #//-------------------------------------------------------------------------------------------#
    #%This is the step function for global level changes e.g. ticks
    #* RUN ONCE PER MODEL PER CALL (I.E. ONCE PER TICK)
    """
    TODO some proper function documentation
    """
    function model_step!(model)
        grid = collect(positions(model))
        tick = model.tick

        saplings = model.saplings
        seedlings = model.seedlings

        sapling_density = model.sapling_density
        seedling_density = model.seedling_density

        n_species = model.n_species

        if model.disturbance_freq > 0
            disturbance_functions.lsp_disturbance(model,
                                                  grid,
                                                  model.disturbance_freq,
                                                  model.disturbed,
                                                  model.max_disturb_size,
                                                  model.nhb_set_ids,
                                                  model.close_nhbs_count,
                                                  model.pcor,
                                                  model.patch_ID,
                                                  n_species,
                                                  seedlings,
                                                  saplings
                                                  )
        end

        if model.external_rain == true
            demog_funcs.external_ldd(model.ext_dispersal_scenario,
                                     grid,
                                     n_species,
                                     model.abundances,
                                     model.external_species,
                                     seedlings)
        end

        for i in 1:length(grid)
            if model.expand[i] == true
                demog_funcs.expand_gap(i, 
                                       model,
                                       grid)
            end

            model.nhb_shade_height[i] = set_get_functions.get_nhb_shade_height(i, 
                                                                               model,
                                                                               grid,
                                                                               range(0, 32, step = 4),
                                                                               model.shell_layers)

            model.nhb_light[i] = set_get_functions.get_light_env(model.nhb_shade_height[i], 
                                                                model.max_heights)

            
        end

        empty_patches = Tuple{Int64, Int64}[]
        for e in empty_positions(model)
            push!(empty_patches, e)
        end 

        empty_patches = Random.shuffle!(empty_patches)
        for p in eachindex(collect(copy(empty_patches)))
            cell_ID =[p]
            demog_funcs.capture_gap(cell_ID, 
                                    model, 
                                    seedlings,
                                    saplings,
                                    model.nhb_light,
                                    model.shade_tolerance,
                                    model.growth_forms,
                                    model.b2_jabowas,
                                    model.b3_jabowas,
                                    model.last_change_tick,
                                    tick,
                                    model.n_changes,
                                    model.pcor)
        end

        count_ep = min(length(collect(empty_positions(model))),
                        length(model.new_agents_list))

        idxs = Random.shuffle(collect(1:count_ep))

        for e in 1:count_ep
            idx = idxs[e]
            pos = empty_patches[idx]
            add_agent!(pos,
                model,
                model.new_agents_list[idx][1], 
                model.new_agents_list[idx][2],
                model.new_agents_list[idx][3], 
                model.new_agents_list[idx][4], 
                model.new_agents_list[idx][5], 
                model.new_agents_list[idx][6], 
                Float64[]
                )
        end

        model.new_agents_list = Any[]

        for p in eachindex(grid)
            a_id = id_in_position(p, model::ABM{<:GridSpaceSingle})
            if a_id !== 0 && model[a_id].patch_here_ID == model.patch_ID[p]
                model[a_id].patch_here_ID = model.patch_ID[p]
            end
        end

        if model.restoration_planting == true && mod(model.tick, model.planting_frequency) == 0
            for i in 1:length(grid)
                saplings[i] .+= model.saplings_to_plant
            end
        end

        seedling_mortality = model.seedling_mortality
        sapling_mortality = model.sapling_mortality
        seedling_transition = model.seedling_transition
        for i in 1:length(grid)
            demog_funcs.regenerate_patch_bank(i,
                                            seedlings, 
                                            saplings,
                                            seedling_mortality,
                                            sapling_mortality,
                                            seedling_transition)
            #TODO may be quicker to have a calc density function
            seedling_density[i] = sum(seedlings[i])
            sapling_density[i] = sum(saplings[i])
        end

        model.max_density = maximum(sapling_density)

        set_get_functions.update_abundances(model, n_species)

        tick += 1
    end
end
