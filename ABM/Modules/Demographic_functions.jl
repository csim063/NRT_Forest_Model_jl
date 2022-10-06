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
        shade_height::Float64,
        height::Float64,
        comp_multiplier::Float64,
        edge_effects::Bool,
        edge_weight::Float64,
        edge_responses::Float64,
        growth_form::Int64,
        dbh::Float64,
        g_jabowas::Float64,
        b2_jabowas::Float64,
        b3_jabowas::Float64,
        max_dbhs::Float64,
        max_heights::Int64
    )
        agent.age += 1

        #*Correct for nhb effects
        competitive_penalty = 1
        if shade_height ≥ height
            comp_val = Float64(abs(Complex(height ./ shade_height) .^ 0.5))
            competitive_penalty = min(1, (comp_multiplier .* comp_val))
        end

        #* Account for edge effects if any
        edge_penalty = 1
        if edge_effects == true
            edge_penalty = (1 - (edge_weight * (1 - edge_responses)))
        end
        
        #* Store growth penalty history for suppression mortality
        growth_reduction = competitive_penalty * edge_penalty
        prepend!(agent.previous_growth, growth_reduction)
        #check if agent list longer than 5 if so select only the first 5
        if length(agent.previous_growth) > 5
            agent.previous_growth = agent.previous_growth[1:5]
        end

        if growth_form == 1
            dbh_increment = ((dbh * g_jabowas * 
                    (1 - (dbh * height) / (max_dbhs * max_heights))) /
                (2.74 + 3 * b2_jabowas * dbh - 4 * b3_jabowas * dbh ^ 2))

            agent.dbh += (dbh_increment * competitive_penalty * edge_penalty)
            agent.height = (1.37 + (b2_jabowas * agent.dbh) - 
                (b3_jabowas * agent.dbh * agent.dbh))

        elseif growth_form == 2
            a_tf_height = 0.05289
            b_tf_height = -0.05695

            hgt_increment = a_tf_height * exp(height * b_tf_height)
            
            agent.height += (hgt_increment * competitive_penalty * edge_penalty)

        else
            error("""The growth form is undefined, 
                     please check species demography data""")
        end
    end

    """ 
    Neighbourhood dispersal
    """
    function nhb_dispersal(
        model,
        seed_production::Int64,
        ldd_disp_frac::Float64,
        r_hgt::Int64,
        DBH::Float64,
        cell_grain::Int64,
        shell_layers::Int64,
        pos::Tuple{Int64,Int64},
        nhbs_id::Vector{Int64},
        shad_h::Vector{Float64},
        seedlings::Vector{Vector{Int64}},
        spec_ID::Int64
    )
        n_seeds = trunc(Int, (rand(Poisson(seed_production)) * 
                    (1 - ldd_disp_frac)))
        
        #? This is a simplified allometric relationship from SORTIE-NZ
        # TODO set constants
        cw = (0.284 .* ((DBH .* 100.0) .^ 0.654))
        shell_width = Int64(min(ceil(cw ./ cell_grain), shell_layers))
        nhb_count = length(collect(nearby_positions(pos, 
                                                model::ABM{<:GridSpaceSingle}, 
                                                shell_width)))
        
        dispersal_nhb_id = nhbs_id[1:nhb_count]


        set_get_functions.assign_seedling(n_seeds::Int64,
                                        dispersal_nhb_id::Vector{Int64},
                                        shad_h::Vector{Float64},
                                        r_hgt::Int64,
                                        seedlings::Vector{Vector{Int64}},
                                        spec_ID::Int64)
    end

    """ 
    Within in patch long distance dispersal
    """
    function ldd_within(
        agent,
        model,
        seed_production::Int64,
        ldd_disp_frac::Float64,
        ldd_dispersal_dist::Int64,
        cell_grain::Int64,
        pos::Tuple{Int64,Int64},
        pcors::Vector{Vector{Tuple{Int64, Int64}}},
        shad_h::Vector{Float64},
        r_hgt::Int64,
        seedlings::Vector{Vector{Int64}},
        spec_ID::Int64
    )
        ldd_seeds = trunc(Int, (rand(Poisson(seed_production)) .* (ldd_disp_frac)))

        #TODO this feels messy and likely slow (SEEMS TO ADD ABOUT 0.4sec per iteration alone)
        for _ in 1:ldd_seeds
            D = trunc(rand(Exponential((ldd_dispersal_dist) ./ cell_grain)));
            
            #TODO change this if you can figure out a way to select just cells at distance not within distance
            # Suggested solution is to alter source code of nearby_positions
            a = collect(nearby_positions(pos, model::ABM{<:GridSpaceSingle}, (D-1)));
            b = collect(nearby_positions(pos, model::ABM{<:GridSpaceSingle}, D));
            D_nhbs = setdiff(b,a);


            target = Int64[]
            if isempty(D_nhbs) .== false
                n = rand(D_nhbs)
                target = [findfirst(isequal([n]), pcors)]
            end

            if length(target) > 0 #seeds that disperse beyond the patch are lost
                target_id = rand(target)
                if shad_h[target_id] <= r_hgt
                    seedlings[target_id][spec_ID] += 1
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
        ext_dispersal_scenario::String,
        grid::Matrix{Tuple{Int64, Int64}},
        n_species::Int64,
        abundances::Vector{Int64},
        external_species::Vector{Float64},
        seedlings::Vector{Vector{Int64}}
    )
        if ext_dispersal_scenario == "equal"
            scalar = trunc((length(grid)./n_species))
            ldd_abundances = fill(scalar, n_species)

        elseif ext_dispersal_scenario == "abundance"
            ldd_abundances = abundances
        else
            error("Unknown external dispersal scenario ($ext_dispersal_scenario), must be 'equal' or 'abundance")
        end
        
        #* Rescale to get the number of seeds to disperse
        crit_min = 0.01 .* length(grid)
        ldd_abundances[ldd_abundances.<crit_min] .= crit_min
        fun = (n::Vector{Int64},p::Vector{Float64})->rand(Binomial(n,p))
        ldd_disperse = map(fun,ldd_abundances, external_species)

        #* Disperse seeds
        #TODO Probably a better way perhaps using the nagents function
        # for s in 1:model.n_species
        for s in eachindex(n_species)
            for i in eachindex(ldd_disperse)#1:length(ldd_disperse)
                for _ in 1:ldd_disperse[i]
                    seedlings[rand(1:length(grid))][s] += 1
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
        cell::Int64,
        ddm::Bool,
        species_ID::Int64,
        base_mortality::Vector{Float64},
        gap_maker::Vector{Int64},
        expand::BitVector,
        previous_species::Vector{Float64},
        previous_height::Vector{Float64},
        a_height::Float64,
        id::Int64,
        age::Float64,
        previous_growth::Vector{Float64},
        supp_tolerance::Vector{Float64},
        supp_mortality::Vector{Float64}
    )
        mort_w = 1

        #* Density dependent mortality
        if ddm == true
            s_nhb = 0
            for i in nearby_ids(agent, model)
                if model[i].species_ID == species_ID
                    s_nhb += 1
                end
            end
            if s_nhb >=4
                mort_w = (1 + ((s_nhb - 4) / 4) * 1)
            end
        end

        #* 1. Baseline-mortality   
        #*   Approximates the standard gap mortality model (see Keane et al. 2001)\
        if rand(Uniform(0.0, 1.0)) < (base_mortality[species_ID] .* mort_w)
            if gap_maker[species_ID] == 1
                expand[cell] = true
            end
            #! NOTE only saving last previous could change to list to record all
            #! but this may impact performance
            previous_species[cell] = species_ID
            previous_height[cell] = a_height

            #model.total_deaths += 1
            kill_agent!(id, model)
            return

        #* 2. Suppression via low growth (p = growth-death if growth < crit-growth [a proportion]).
        #* trees older than 10 years have a chance of competition based death
        elseif (
            age > 10.0 && 
            mean(previous_growth) < supp_tolerance[species_ID] &&
            rand(Uniform(0.0, 1.0)) < supp_mortality[species_ID]
            )
            if gap_maker[species_ID] == 1
                expand[cell] = true
            end
            #! NOTE only saving last previous could change to list to record all
            #! but this may impact performance
            previous_species[cell] = species_ID
            previous_height[cell] = a_height

            #model.total_deaths += 1
            kill_agent!(id, model)
            return
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
                cell = model[i].patch_here_ID
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
        cell_ID::Vector{Int64},
        model,
        seedlings::Vector{Vector{Int64}},
        saplings::Vector{Vector{Int64}},
        nhb_light::Vector{Float64},
        shade_tolerance::Vector{Float64},
        growth_forms::Vector{Int64},
        b2_jabowas::Vector{Float64},
        b3_jabowas::Vector{Float64},
        last_change_tick::Vector{Int64},
        tick::Int64,
        n_changes::Vector{Int64},
        pcor::Vector{Vector{Tuple{Int64, Int64}}}
    )
        seedlings = seedlings[cell_ID]
        saplings = saplings[cell_ID]    

        #? 0.25 is probability of one sapling becoming an adult
        if sum(saplings[1]) >= 1 && rand(Uniform(0, 1)) < (1.0 - (0.25 ^ sum(saplings[1])))
            weights = floor.(Int, (abs.(nhb_light[cell_ID] .- 
                                        shade_tolerance)) * 100)

            regenbank_wgt = saplings[1] .* weights

            new_species_id = distribution_functions.lottery(regenbank_wgt, true)
            new_species_id = new_species_id !== nothing ? new_species_id : rand(1:length(seedlings))
            if growth_forms[new_species_id] == 1
                #? dbh (in m) of 0.01 m (1 cm) + noise (0, 0.01)
                dbh = 0.01 .+ rand(Uniform(0, 0.01))
                #TODO CONFIRM IN NETLOGO IT USES ONLY b2-jabowa I THINK THIS IS A MISTAKE 
                #TODO HERE I USE b3 AS IS USED IN INITIALISATION
                b2 = b2_jabowas[new_species_id]
                b3 = b3_jabowas[new_species_id]
                height = 1.37 .+ (b2 .* dbh) .- (b3 .* dbh .* dbh)
                age = 1.0
            elseif growth_forms[new_species_id] == 2
                #TODO CONFIRM DBH OF GROWTH FORM 2 IN NETLOGO THEY ARE NOT DEFINED
                dbh = 0.01 .+ rand(Uniform(0, 0.01)) #!LIKELY WRONG
                height = 1.5 .+ rand(Uniform(0, 0.1))
                age = 1.0
            else
                error("The growth form for $new_species_id is undefined, please check species demography data")
            end

            last_change_tick[cell_ID] .= tick
            n_changes[cell_ID] .+= 1

            #* Add new adult tree as agent
            new_agent = [new_species_id,
                        cell_ID[1],
                        growth_forms[new_species_id],
                        height,
                        dbh,
                        age]

            push!(model.new_agents_list, new_agent)

            #* Empty the regeneration bank
            seedlings = seedlings .- seedlings
            saplings = saplings .- saplings
        end
    end

    """
    Function to adjust seedling and sapling values in each patch
    """
    function regenerate_patch_bank(
        cell_ID::Int64,
        seedlings::Vector{Vector{Int64}},
        saplings::Vector{Vector{Int64}},
        seedling_mortality::Vector{Float64},
        sapling_mortality::Vector{Float64},
        seedling_transition::Vector{Float64}
    )
        dead_saplings = rand.(Binomial.(abs.(saplings[cell_ID]), sapling_mortality))
        dead_seedlings = rand.(Binomial.(abs.(seedlings[cell_ID]), seedling_mortality))

        saplings[cell_ID] = saplings[cell_ID] - dead_saplings
        saplings[cell_ID] = replace(x-> x < 0 ? 0 : x, saplings[cell_ID])
        seedlings[cell_ID] = seedlings[cell_ID] - dead_seedlings
        seedlings[cell_ID] = replace(x-> x < 0 ? 0 : x, seedlings[cell_ID])

        #* New saplings
        new_saplings = rand.(Binomial.(abs.(seedlings[cell_ID]), seedling_transition))

        saplings[cell_ID] = saplings[cell_ID] + new_saplings
        seedlings[cell_ID] = seedlings[cell_ID] + new_saplings
    end
end