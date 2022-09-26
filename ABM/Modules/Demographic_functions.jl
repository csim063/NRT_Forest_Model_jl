## TODO Documentation

module demog_funcs
    using Agents
    using StatsBase
    using Distributions

    include("Helper_functions.jl")

    """
    Grow each tree following basic equations of Botkin et al. (1972) JABOWA
    down-weight based on relative hgt of cell (following Dislich et al. 2009)
    """
    function grow(
        agent,
        model
    )
        spec_num = agent.species_ID
        agent.age += 1

        #*Correct for nhb effects
        competitive_penalty = 1
        if model.nhb_shade_height[spec_num] ≥ agent.height
            comp_val = Float64(abs(Complex(agent.height / model.nhb_shade_height[spec_num]) ^ 0.5))
            competitive_penalty = min(1, (model.comp_multiplier * comp_val))
        end

        #* Account for edge effects if any
        edge_penalty = 1
        if model.edge_effects == true
            edge_penalty = (1 - (model.edge_weight[spec_num] * (1 - model.edge_responses[spec_num])))
        end
        
        #* Store growth penalty history for suppression mortality
        growth_reduction = competitive_penalty * edge_penalty
        prepend!(agent.previous_growth, growth_reduction)
        #check if agent list longer than 5 if so select only the first 5
        if length(agent.previous_growth) > 5
            agent.previous_growth = agent.previous_growth[1:5]
        end

        if agent.growth_form == 1
            dbh_increment = ((agent.dbh * model.g_jabowas[spec_num] * 
                    (1 - (agent.dbh * agent.height) / (model.max_dbhs[spec_num] * model.max_heights[spec_num]))) /
                (2.74 + 3 * model.b2_jabowas[spec_num] * agent.dbh - 4 * model.b3_jabowas[spec_num] * agent.dbh ^ 2))

            agent.dbh += (dbh_increment * competitive_penalty * edge_penalty)
            agent.height = (1.37 + (model.b2_jabowas[spec_num] * agent.dbh) - 
                (model.b3_jabowas[spec_num] * agent.dbh * agent.dbh))

        elseif agent.growth_form == 2
            a_tf_height = 0.05289
            b_tf_height = -0.05695

            hgt_increment = a_tf_height * exp(agent.height * b_tf_height)
            
            agent.height += (hgt_increment * competitive_penalty * edge_penalty)

        else
            error("""The growth form $agent.growth_form is undefined, 
                     please check species demography data""")
        end
    end

    """ 
    Neighbourhood dispersal
    """
    function nhb_dispersal(
        agent,
        model
    )
        if model.seed_prod[agent.species_ID] > 0
            n_seeds = trunc(Int, (rand(Poisson(model.seed_prod[agent.species_ID])) * 
                        (1 - model.ldd_dispersal_fracs[agent.species_ID])))
            
            r_hgt = model.regen_heights[agent.species_ID]
            
            #? This is a simplified allometric relationship from SORTIE-NZ
            cw = (0.284 * ((agent.dbh * 100) ^ 0.684))
            shell_width = Int(min(ceil(cw / model.cell_grain), length(model.nhb_set[agent.species_ID])))
            dispersal_nhb = model.nhb_set[agent.species_ID][1:shell_width] #! SEEMS WRONG AS SPECIES ID SEEMS ARBITRARY AS A SELECTOR

            for _ in 1:n_seeds
                rand_cell = rand(rand(dispersal_nhb)) #THIS IS A (XCOR, YCOR)
                #TODO try replace findfirst
                rand_cell_ID = findfirst(x->x==[rand_cell], model.pcor)[1]
                
                if model.nhb_shade_height[rand_cell_ID] ≤ r_hgt
                    model.seedlings[rand_cell_ID][agent.species_ID] += 1
                end
            end
        end
    end

    """ 
    Within in patch long distance dispersal
    """
    function ldd_within(
        agent,
        model
    )
        ldd_seeds = trunc(Int, (rand(Poisson(model.seed_prod[agent.species_ID])) * 
                    (model.ldd_dispersal_fracs[agent.species_ID])))

        #TODO this feels messy and likely slow (SEEMS TO ADD ABOUT 0.4sec per iteration alone)
        for _ in 1:ldd_seeds
            D = trunc(rand(Exponential((model.ldd_dispersal_dist[agent.species_ID]) / model.cell_grain)));
            
            #TODO change this if you can find a way to select just cells at distance not within distance
            # MAYBE x[findall(y -> y == 10, x)[1]] could be useful
            # Suggested solution is to alter source code of nearby_positions
            a = collect(nearby_positions(agent.pos, model::ABM{<:GridSpaceSingle}, (D-1)));
            b = collect(nearby_positions(agent.pos, model::ABM{<:GridSpaceSingle}, D));
            D_nhbs = setdiff(b,a);

            target = isempty(D_nhbs) ? Int64[] : findall(x->x==[rand(D_nhbs)], model.pcor)
            if length(target) > 0 #seeds that disperse beyond the patch are lost
                target_id = rand(target)
                if model.nhb_shade_height[target_id] <= model.regen_heights[agent.species_ID]
                    model.seedlings[target_id][agent.species_ID] += 1
                end
            end
        end
    end

    """
    Long-distance dispersal from *beyond* the plot
    Compute total LDD seed rain and then disperse it across patches
    NOTE THIS IS ONLY RUN ONCE PER TICK NOT ONCE PER AGENT PER TICK
    """
    function external_ldd(
        model,
        ext_dispersal_scenario,
        grid
    )
        if ext_dispersal_scenario == "equal"
            scalar = trunc((length((positions(model)))/model.n_species))
            ldd_abundances = fill(scalar, model.n_species)

        elseif ext_dispersal_scenario == "abundance"
            ldd_abundances = model.abundances
        else
            error("Unknown external dispersal scenario ($ext_dispersal_scenario), must be 'equal' or 'abundance")
        end
        
        #* Rescale to get the number of seeds to disperse
        crit_min = 0.01 * length((positions(model)))
        ldd_abundances[ldd_abundances.<crit_min] .= crit_min
        fun = (n,p)->rand(Binomial(n,p))
        ldd_disperse = map(fun,ldd_abundances, model.external_species)

        #* Disperse seeds
        #TODO Probably a better way perhaps using the nagents function
        for s in 1:model.n_species
            for i in eachindex(ldd_disperse)#1:length(ldd_disperse)
                for _ in 1:ldd_disperse[i]
                    model.seedlings[rand(1:length(grid))][s] += 1
                end
            end
        end
    end

    """
    Function to determine impact of herbivory
    """
    function herbivore_effect(
        agent,
        model
    )
        cell = agent.patch_here_ID

        model.seedlings[cell] = round.(Int64, model.seedlings[cell] .* model.herbivory_amount)
        if model.saplings_eaten == true
            model.saplings[cell] = round.(Int64, model.saplings[cell] .* model.herbivory_amount)
        end 
    end

    """
    Function to add species-specific local effects on the seedling bank
    """
    function thin_regenbank(
        agent,
        model
    )
        inhibit = model.seedling_inhibition[agent.species_ID] 
        if inhibit > 0 && inhibit < 1
            #cell = findfirst(x->x==[agent.pos], model.pcor)[1]
            cell = agent.patch_here_ID
            model.seedlings[cell] = round.(Int64, model.seedlings[cell] .* (1 - inhibit))
        end
    end

    """
    Function to get impact of tree-fern litter fall on saplings
    """
    function macro_litter_fall(
        agent,
        model
    )
        if rand(Uniform(0, 1)) < model.macro_litter_effect
            # cell = findfirst(x->x==[agent.pos], model.pcor)[1] 
            cell = agent.patch_here_ID
            spp_die = rand(1:model.n_species)
            model.saplings[cell][spp_die] -= 1
        end
    end

    """
    Kill trees at background mortality rate
    """
    function death(
        agent,
        model, 
    )
        mort_w = 1
        # cell = findfirst(x->x==[agent.pos], model.pcor)[1] 
        cell = agent.patch_here_ID

        #* Density dependent mortality
        if model.ddm == true
            s_nhb = 0
            for i in nearby_ids(agent, model)
                if model[i].species_ID == agent.species_ID
                    s_nhb += 1
                end
            end
            if s_nhb >=4
                mort_w = (1 + ((s_nhb - 4) / 4) * 1)
            end
        end

        #* 1. Baseline-mortality   
        #*   Approximates the standard gap mortality model (see Keane et al. 2001)\
        if rand(Uniform(0, 1)) < (model.base_mortality[agent.species_ID] * mort_w)
            if model.gap_maker[agent.species_ID] == 1
                model.expand[cell] = true
            end
            #! NOTE only saving last previous could change to list to record all
            #! but this may impact performance
            model.previous_species[cell] = agent.species_ID
            model.previous_height[cell] = agent.height

            kill_agent!(agent.id, model)

        #* 2. Suppression via low growth (p = growth-death if growth < crit-growth [a proportion]).
        #* trees older than 10 years have a chance of competition based death
        elseif (
            agent.age > 10.0 && 
            mean(agent.previous_growth) < model.supp_tolerance[agent.species_ID] &&
            rand(Uniform(0, 1)) < model.supp_mortality[agent.species_ID]
            )
            if model.gap_maker[agent.species_ID] == 1
                model.expand[cell] = true
            end
            #! NOTE only saving last previous could change to list to record all
            #! but this may impact performance
            model.previous_species[cell] = agent.species_ID
            model.previous_height[cell] = agent.height

            kill_agent!(agent.id, model)
        end
    end

    """
    Grow the gap if the dying species is a 'gap maker': a cone of neighbouring cells is affected
    """
    #! NOTE THIS FUNCTION IS DIFFERENT TO NETLOGO AS IT USES PREVIOUS HEIGHT AS AGENT IS DEAD
    #! AT TIME OF CALLING (NETLOGO ALWAYS HAS A HEIGHT OF 0) 
    #! AND A CONE ANGLE OF ZERO SO JUST SELECTS CURRENT PATCH NO GROWTH
    function expand_gap(
        cell_ID,
        model,
        grid
    )
        h = trunc(Int64, (model.previous_height[cell_ID] / model.cell_grain) + 1)

        nearby_trees = nearby_ids(grid[cell_ID], model::ABM{<:GridSpaceSingle}, h)
        ax = rand(1:2)
        fall_dir = model.pcor[cell_ID][1][ax]
        for i in nearby_trees 
            if model[i].pos[ax] == fall_dir
                cell = findfirst(x->x==[model[i].pos], model.pcor)[1] 
                model.previous_species[cell] = model[i].species_ID
                model.previous_height[cell] = model[i].height

            kill_agent!(model[i].id, model)
            end
        end
    end

    """
    Have a new tree grow from sapling in a gap
    """
    function capture_gap(
        cell_ID,
        model
    )
        saplings = model.saplings[cell_ID]    

        #? 0.25 is probability of one sapling becoming an adult
        if sum(saplings)[1] >= 1 && rand(Uniform(0, 1)) < (1.0 - (0.25 ^ sum(saplings)[1]))
            weights = floor.(Int, (abs.(model.nhb_light[cell_ID] .- 
                                        model.shade_tolerance)) * 100)

            regenbank_wgt = saplings .* weights

            new_species_id = distribution_functions.lottery(regenbank_wgt, true)
            if model.growth_forms[new_species_id] == 1
                #? dbh (in m) of 0.01 m (1 cm) + noise (0, 0.01)
                dbh = 0.01 + rand(Uniform(0, 0.01))
                #TODO CONFIRM IN NETLOGO IT USES ONLY b2-jabowa I THINK THIS IS A MISTAKE 
                #TODO HERE I USE b3 AS IS USED IN INITIALISATION
                b2 = model.b2_jabowas[new_species_id]
                b3 = model.b3_jabowas[new_species_id]
                height = 1.37 + (b2 * dbh) - (b3 * dbh * dbh)
                age = 1.0
            elseif model.growth_forms[new_species_id] == 2
                #TODO CONFIRM DBH OF GROWTH FORM 2 IN NETLOGO THEY ARE NOT DEFINED
                dbh = 0.01 + rand(Uniform(0, 0.01)) #!LIKELY WRONG
                height = 1.5 + rand(Uniform(0, 0.1))
                age = 1.0
            else
                error("The growth form for $new_species_id is undefined, please check species demography data")
            end

            #TODO could easily change last change tick into list of all ticks where changes occur
            model.last_change_tick[cell_ID] .= model.tick
            model.n_changes[cell_ID] .+= 1

            #* Add new adult tree as agent
            #! ABSOLUTELY NO CLUE WHY I NEED TO CREATE A NEW TUPLE AND DIVE 3 LAYERS IN TO EACH VARIABLE
            posit = (model.pcor[cell_ID][1][1][1], model.pcor[cell_ID][1][1][2])

            add_agent!(posit, model, 
                new_species_id, 
                cell_ID[1],
                model.growth_forms[new_species_id], 
                height, 
                dbh, 
                age, 
                Float64[]
                )

            #* Empty the regeneration bank
            model.seedlings[cell_ID] = model.seedlings[cell_ID] - model.seedlings[cell_ID]
            model.saplings[cell_ID] = model.saplings[cell_ID] - model.saplings[cell_ID]
        end
    end

    """
    Function to adjust seedling and sapling values in each patch
    """
    function regenerate_patch_bank(
        cell_ID,
        model
    )
        dead_saplings = rand.(Binomial.(abs.(model.saplings[cell_ID]), model.sapling_mortality))
        dead_seedlings = rand.(Binomial.(abs.(model.seedlings[cell_ID]), model.seedling_mortality))

        model.saplings[cell_ID] = model.saplings[cell_ID] - dead_saplings
        model.saplings[cell_ID] = replace(x-> x < 0 ? 0 : x, model.saplings[cell_ID])
        model.seedlings[cell_ID] = model.seedlings[cell_ID] - dead_seedlings
        model.seedlings[cell_ID] = replace(x-> x < 0 ? 0 : x, model.seedlings[cell_ID])

        #* New saplings
        new_saplings = rand.(Binomial.(abs.(model.seedlings[cell_ID]), model.seedling_transition))

        model.saplings[cell_ID] = model.saplings[cell_ID] + new_saplings
        model.seedlings[cell_ID] = model.seedlings[cell_ID] + new_saplings
    end
end