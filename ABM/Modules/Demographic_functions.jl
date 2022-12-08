"""
Module contains the core functions which enable agents to grow, reproduce, expand their ranges, and
die.
"""

module demog_funcs
    using Agents
    using StatsBase
    using Distributions

    #* CUSTOM MODULES
    include("Helper_functions.jl")

    #//-------------------------------------------------------------------------------------------#
    #% GROW
    """
    # Grow
    Have each tree grow, i.e. increase their age, height and diameter at breast height. Every tree
    increases their age by 1 every tick and their height and DBH are incremented based upon their
    growth form (either tree or tree-fern) and basic demography equations from Botkin et *al* (1972)
    and Brock et *al* (2019). The incremental growth in height and DBH is down-weighted based on the
    relative amount of shading (due to the height) from neighbouring cells (following Dislich et 
    *al* 2009).
    ## References:
    - Botkin, D.B., Janak, J.F., Wallis, J.R., 1972. Some ecological consequences of acomputer model 
    of forest growth. J. Ecol. 60, 849-872.
    - Dislich, C., Günter, S., Homeier, J., Schröder, B., & Huth, A. (2009). Simulating Forest 
    Dynamics of a Tropical Montane Forest in South Ecuador. Erdkunde, 63(4), 347-364.
    - Brock, J. M. R., Morales, N. S., Burns, B. R., & Perry, G. L. W. (2019). The hare, tortoise 
    and crocodile revisited: Tree fern facilitation of conifer persistence and angiosperm growth in 
    simulated forests. The Journal of Ecology, 108(3), 969-981.
    ## Arguments:
    - `agent`: Agent object of type `GridAgent` with properties matching those defined by the 
    current model. Note agent must exist (i.e. not have died).
    - `shade_height::Float64`: Mean weighted height of neighbouring trees as calculated by
    `get_nhb_shade_height()`.
    - `height::Float64`: Height in meters of tree.
    - `comp_multiplier::Float64`: Constant impacting the strength of competition that each species 
    experiences. This constant is α in α * (height / shade_height) ^ 0.5 based on Dislich et 
    *al* 2009.
    - `edge_effects::Bool`: Boolean determining whether trees are impacted by edge effects or not.
    - `edge_weight::Float64`: The strength of the edge effect relative to the distance from the 
    edge.
    - `edge_responses::Float64`: Value between 0 and 1 representing how well a species is adapted to 
    edge environments
    - `growth_form::Int64`: Growth type of agent, 1 = trees; 2 = tree ferns.
    - `dbh::Float64`: Diameter at breast height in meters of tree.
    - `g_jabowas::Float64`: Species-specific allometric constant defined by Botkin et *al* (1972).
    - `b2_jabowas::Float64`: Species-specific allometric constant defined by Botkin et *al* (1972)
    defined as `2 * (max_height - 1.37) / max_dbh` where 1.37 is the smallest adult tree height.
    - `b3_jabowas::Float64`: Species-specific allometric constant defined by Botkin et *al* (1972)
    defined as `(max_height - 1.37) / max_dbh ^ 2` where 1.37 is the smallest adult tree height.
    - `max_dbhs::Float64`: Species-specific maximum attainable diameter at breast height (m).
    - `max_heights::Int64`: Species-specific maximum attainable height (m).
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

        #% CALCULATE COMPETITIVE PENALTY FROM NEIGHBOURS------------#
        competitive_penalty = 1
        if shade_height ≥ height
            comp_val = Float64(abs(Complex(height ./ shade_height) .^ 0.5))
            competitive_penalty = min(1, (comp_multiplier .* comp_val))
        end

        #% CALCULATE EDGE PENALTY-----------------------------------#
        edge_penalty = 1
        if edge_effects == true
            edge_penalty = (1 - (edge_weight * (1 - edge_responses)))
        end
        
        #% CALCULATE OVERALL GROWTH PENALTY-------------------------#
        #* Store growth penalty history for use
        #* later in suppression mortality
        growth_reduction = competitive_penalty * edge_penalty
        prepend!(agent.previous_growth, growth_reduction)

        #* Only the last 5 ticks worth of growth suppression are 
        #* stored
        if length(agent.previous_growth) > 5
            agent.previous_growth = agent.previous_growth[1:5]
        end

        #% INCREMENT TREE HEIGHT AND DBH----------------------------#
        # For trees
        if growth_form == 1
            #? Equations are based on Botkin et al (1972)
            dbh_increment = ((dbh * g_jabowas * 
                    (1 - (dbh * height) / (max_dbhs * max_heights))) /
                (2.74 + 3 * b2_jabowas * dbh - 4 * b3_jabowas * dbh ^ 2))

            #! Note change of scoping when assigning new values. This is done
            #! to ensure values are updated in the global pool
            agent.dbh += (dbh_increment * competitive_penalty * edge_penalty)
            agent.height = (1.37 + (b2_jabowas * agent.dbh) - 
                (b3_jabowas * agent.dbh * agent.dbh))
        
        # For tree ferns
        elseif growth_form == 2
            #? Equations and values are based on Brock et al (2019)
            a_tf_height = 0.05289
            b_tf_height = -0.05695

            hgt_increment = a_tf_height * exp(height * b_tf_height)
            
            agent.height += (hgt_increment * competitive_penalty * edge_penalty)

        else
            error("""The growth form is undefined, 
                     please check species demography data""")
        end
    end

    #//-------------------------------------------------------------------------------------------#
    #% NEIGHBOURHOOD DISPERSAL
    """ 
    # Neighbourhood dispersal
    Seeds, represented as seedlings, are dispersed from target tree. Seeds are assumed to land 
    under the crown of the target tree. The crown size of the tree is computed based about its 
    diameter at breast height following a simplified allometric relationship from SORTIE-NZ
    (Kunstler et *al*., 2009).
    ## References:
    - Kunstler, G., Coomes, D.A., Canham, C.D., 2009. Size-dependence of growth and mortality 
    influence the shade tolerance of trees in a lowland temperate rain forest. J. Ecol. 97, 685-695. 
    ## Arguments:
    - `model`: The AgentBasedModel object defining the current model. This object is usually 
    created using `Agents.ABM()`.
    - `seed_production::Int64`: Species-specific constant value used within a Poisson distribution
    to determine the number of seeds produced by an agent.
    - `ldd_disp_frac::Float64`: Proportion of seedlings that undergo long-distance dispersal and are
    therefore unavailable for neighbourhood dispersal.
    - `r_hgt::Int64`: Regeneration height in meters for selected species.
    - `dbh::Float64`: Diameter at breast height in meters of tree.
    - `cell_grain::Int64`: The size of an individual cell in the model (cell_grain x cell_grain).
    - `shell_layers::Int64`: Maximum distance at which neighbours may be considered to possibly fall
    under the crown.
    - `pos::Tuple{Int64,Int64}`: Coordinates of tree/agent.
    - `nhbs_id::Vector{Int64}`: Patch IDs for all neighbouring cells of tree/agent.
    - `shad_h::Vector{Float64}`: Mean shade heights as calculated by `get_nhb_shade_height()` for 
    all cells in model grid.
    - `seedlings::Vector{Vector{Int64}}`: Number of seedlings for each species for every cell in the
    model grid.
    - `spec_ID::Int64`: Selected species ID.
    """
    function nhb_dispersal(
        model,
        seed_production::Int64,
        ldd_disp_frac::Float64,
        r_hgt::Int64,
        dbh::Float64,
        cell_grain::Int64,
        shell_layers::Int64,
        pos::Tuple{Int64,Int64},
        nhbs_id::Vector{Int64},
        shad_h::Vector{Float64},
        seedlings::Vector{Vector{Int64}},
        spec_ID::Int64
    )
        #% GENERATE SEEDS TO BE DISPERSED---------------------------#
        n_seeds = trunc(Int, (rand(Poisson(seed_production)) * 
                    (1 - ldd_disp_frac)))
        
        #% CALCULATE CROWN WIDTH AND NEIGHBOURS---------------------#
        #? This is a simplified allometric relationship 
        #?from SORTIE-NZ
        cw = (0.284 .* ((dbh .* 100.0) .^ 0.654))
        shell_width = Int64(min(ceil(cw ./ cell_grain), shell_layers))
        nhb_count = length(collect(nearby_positions(pos, 
                                                model::ABM{<:GridSpaceSingle}, 
                                                shell_width)))
        
        #% ASSIGN LOCALLY DISPERSED SEEDS---------------------------#                                        
        dispersal_nhb_id = nhbs_id[1:nhb_count]
        set_get_functions.assign_seedling(n_seeds::Int64,
                                        dispersal_nhb_id::Vector{Int64},
                                        shad_h::Vector{Float64},
                                        r_hgt::Int64,
                                        seedlings::Vector{Vector{Int64}},
                                        spec_ID::Int64)
    end

    #//-------------------------------------------------------------------------------------------#
    #% WITHIN PATCH LONG DISTANCE DISPERSAL
    """ 
    # Within patch long distance dispersal
    Disperse a fraction of seeds beyond the parent crown but within the same habitat fragment. The
    disperesed seeds are added to the seedlings bank for the appropriate grid cell.
    ## Arguments:
    - `model`: The AgentBasedModel object defining the current model. This object is usually 
    created using `Agents.ABM()`.
    - `seed_production::Int64`: Species-specific constant value used within a Poisson distribution
    to determine the number of seeds produced by an agent.
    - `ldd_disp_frac::Float64`: Proportion of seedlings that undergo long-distance dispersal and are
    therefore unavailable for neighbourhood dispersal.
    - `ldd_dispersal_dist::Int64`: Species-specific mean dispersal distance
    - `cell_grain::Int64`: The size of an individual cell in the model (cell_grain x cell_grain).
    - `pos::Tuple{Int64,Int64}`: Coordinates of tree/agent.
    - `pcors::Vector{Vector{Tuple{Int64, Int64}}}`: List of coordinates for all cells in model grid.
    - `shad_h::Vector{Float64}`: Mean shade heights as calculated by `get_nhb_shade_height()` for 
    all cells in model grid.
    - `r_hgt::Int64`: Regeneration height in meters for selected species.
    - `seedlings::Vector{Vector{Int64}}`: Number of seedlings for each species for every cell in the
    model grid.
    - `spec_ID::Int64`: Selected species ID.
    """
    function ldd_within(
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
        #% GENERATE SEEDS TO BE DISPERSED---------------------------#
        ldd_seeds = trunc(Int, (rand(Poisson(seed_production)) 
                                            .* (ldd_disp_frac)))

        #% LOOP THROUGH ALL SEEDS AND ASSIGN TO A CELL--------------#
        #Threads.@threads for _ in 1:ldd_seeds
        for _ in 1:ldd_seeds
            #Calculate appropriate distance
            D = trunc(Int, rand(Exponential((ldd_dispersal_dist) ./ cell_grain)));
            
            #TODO Move this documentation to the top of the function
            #* Find cells at distance D only. To do this currently we find all the
            #* cells within distance D and remove the ones also within distance D-1
            D_nhbs = set_get_functions.positions_at(pos::Tuple{Int64,Int64}, 
                                                    model::ABM{<:GridSpaceSingle}, 
                                                    D::Int64,)

            ## If there are avaible cells select one to be the target
            target = Int64[]
            if isempty(D_nhbs) .== false
                n = rand(D_nhbs)
                target = [findfirst(isequal([n]), pcors)]
            end

            if length(target) > 0
                target_id = rand(target)
                #* If regeneration height is larger than targets
                #* mean shading height assign the seed
                if shad_h[target_id] <= r_hgt
                    seedlings[target_id][spec_ID] += 1
                end
            end
        end
    end

    #//-------------------------------------------------------------------------------------------#
    #% EXTERNAL TO PATCH LONG DISTANCE DISPERSAL
    """
    # External long distance dispersal
    Long-distance dispersal (or seedling rain) from beyond the modelled patch. Each grid cell has a 
    chance of receiving a seedling from outside the simulated fragment. The probability of 
    receiving a seedling is the larger of the species proportional abundance in the grid or a 
    minimum chance of 1%
    ## Arguments:
    - `ext_dispersal_scenario::String`: Whether seeds are equally divided between species 
    (`"equal"`) or based upon current abundances (`"abundance"`).
    - `grid::Matrix{Tuple{Int64, Int64}}`: Matrix containing all of the coordinates of cells in the
    current model.
    - `n_species::Int64`: Number of total species.
    - `abundances::Vector{Int64}`: Count of total trees of each species in the model.
    - `external_species::Vector{Float64}`: Species-specific probability of a seed being actively
    dispersed from beyond the modelled habitat patch.
    - `seedlings::Vector{Vector{Int64}}`: Number of seedlings for each species for every cell in the
    model grid.
    """
    function external_ldd(
        ext_dispersal_scenario::String,
        grid::Matrix{Tuple{Int64, Int64}},
        n_species::Int64,
        abundances::Vector{Int64},
        external_species::Vector{Float64},
        seedlings::Vector{Vector{Int64}}
    )
        #% CALCULATE ABUNDANCE--------------------------------------#
        #* Abundance is defined based on the dispersal scenario selected
        if ext_dispersal_scenario == "equal"
            scalar = trunc((length(grid)./n_species))
            ldd_abundances = fill(scalar, n_species)
        elseif ext_dispersal_scenario == "abundance"
            ldd_abundances = abundances
        else
            error("Unknown external dispersal scenario ($ext_dispersal_scenario), 
            must be 'equal' or 'abundance")
        end
        
        #% RESCALE AND CALCULATE NUMBER OF SEEDS--------------------#
        crit_min = 0.01 .* length(grid)
        ldd_abundances[ldd_abundances.<crit_min] .= crit_min
        #* Anonymous function determines number of successful seed
        #* dispersals based on number of trails and dispersal probability 
        fun = (n::Vector{Int64},p::Vector{Float64})->rand(Binomial(n,p))
        ldd_disperse = map(fun,ldd_abundances, external_species)

        #% DISPERSE SEEDS-------------------------------------------#
        Threads.@threads for s in eachindex(n_species)
            for i in eachindex(ldd_disperse)
                for _ in 1:ldd_disperse[i]
                    seedlings[rand(1:length(grid))][s] += 1
                end
            end
        end
    end

    #//-------------------------------------------------------------------------------------------#
    #% HERBIVORY
    """
    # Calculate herbivore effect
    Reduce the seedling and sapling banks by the fraction that suffer herbivory. Herbivory is 
    explicitly considered as component of juvenile mortality, hence why only seedlings and saplings 
    are impacted.
    ## Arguments:
    - `agent`: Agent object of type `GridAgent` with properties matching those defined by the 
    current model. Note agent must exist (i.e. not have died).
    - `model`: The AgentBasedModel object defining the current model. This object is usually 
    created using `Agents.ABM()`.
    """
    function herbivore_effect(
        agent,
        model
    )
        cell = agent.patch_here_ID

        ## Reduce seedlings
        model.seedlings[cell] = round.(Int64, model.seedlings[cell] .* model.herbivory_amount)

        ## If being considered reduce saplings
        if model.saplings_eaten == true
            model.saplings[cell] = round.(Int64, model.saplings[cell] .* model.herbivory_amount)
        end 
    end

    #//-------------------------------------------------------------------------------------------#
    #% SEEDLING BANK EFFECTS
    """
    # Inhibit seedling bank
    Implements species-specific local effects reducing the seedlings in the target patch.
    ## Arguments:
    - `agent`: Agent object of type `GridAgent` with properties matching those defined by the 
    current model. Note agent must exist (i.e. not have died).
    - `model`: The AgentBasedModel object defining the current model. This object is usually 
    created using `Agents.ABM()`.
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

    #//-------------------------------------------------------------------------------------------#
    #% LITTER FALL
    """
    # Macro-litterfall damage
    Calculate the damage on saplings resulting from macro-litterfall (primarily under tree ferns).
    This function simply tests whether a macro-litterfall event has occured based on a probability
    and if true removes a single sapling from one a random species in the grid cell
    ## Arguments:
    - `agent`: Agent object of type `GridAgent` with properties matching those defined by the 
    current model. Note agent must exist (i.e. not have died).
    - `model`: The AgentBasedModel object defining the current model. This object is usually 
    created using `Agents.ABM()`.
    """
    function macro_litter_fall(
        agent,
        model
    )
        if rand(Uniform(0, 1)) < model.macro_litter_effect
            cell = agent.patch_here_ID
            spp_die = rand(1:model.n_species)
            if model.saplings[cell][spp_die] ≥ 1
                model.saplings[cell][spp_die] -= 1
            end
        end
    end

    #//-------------------------------------------------------------------------------------------#
    #% DEATH
    """
    # Background mortality
    Function applies a baseline mortality to all agents with a given proportion, as well as
    accounting for age and competitive suppression and density dependent mortality. The function 
    is broken into two broad sections. 1) Simple baseline mortality which is simply a probability
    of a mortality event happening in a tick weighted by debnsity dependent mortality if this is 
    selected. 2) Suppression via low growth and aging mortality.
    ## Arguments:
    - `agent`: Agent object of type `GridAgent` with properties matching those defined by the 
    current model. Note agent must exist (i.e. not have died).
    - `model`: The AgentBasedModel object defining the current model. This object is usually 
    created using `Agents.ABM()`.
    - `cell::Int64`: Patch ID of cell where target tree is located
    - `ddm::Bool`: Boolean indicating whether density dependent mortality should be accounted for 
    or not
    - `species_ID::Int64`: Selected species ID.
    - `base_mortality::Vector{Float64}`: Species specific background mortaility rate. Provided as a
    probability (0-1) of a tree dying in anyone tick.
    - `gap_maker::Int64`: Species specific property indicating whether a species is capable 
    of creating a forest gap (value = 1) or not (value = 0).
    - `expand::BitVector`: Flag indicating that a patch should be checked for gap expansion 
    (value = 1) or not (value = 0). See `expand_gap()`.
    - `previous_height::Vector{Float64}`: Final height of the last tree to have occupied every cell 
    in the current model grid.
    - `a_height::Float64`: Height in meters of the tree.
    - `id::Int64`: Agent ID of the selected tree/agent
    - `age::Float64`: Age (number of ticks alive) of tree
    - `previous_growth::Vector{Float64}`: Growth penalty list for the last 5 ticks worth of growth 
    suppression experienced by the selected tree
    - `supp_tolerance::Vector{Float64}`: Species specific suppression tolerance
    - `supp_mortality::Vector{Float64}`: Species specific probability of mortality due to 
    suppression 
    - `grass::Bool`: Boolean indicating whether grass is present in the model or not
    - `grass_invasion_prob::Float64`: Probability of grass invading a patch
    - `pests::Bool`: Boolean indicating whether pests are present in the model or not
    - `pest_mortality::Float64`: Probability of mortality due to pests (this is the mean in a
        normal distribution)
    - `pest_var::Float64`: Variance of the normal distribution used to calculate the probability
    """
    function death(
        agent,
        model, 
        cell::Int64,
        ddm::Bool,
        species_ID::Int64,
        base_mortality::Vector{Float64},
        gap_maker::Int64,
        expand::BitVector,
        previous_height::Vector{Float64},
        a_height::Float64,
        id::Int64,
        age::Float64,
        previous_growth::Vector{Float64},
        supp_tolerance::Vector{Float64},
        supp_mortality::Vector{Float64},
        grass::Bool,
        grass_invasion_prob::Float64,
        pests::Bool,
        pest_mortality::Float64,
        pest_var::Float64,
        weather_adjustments::Float64
    )
        ## Define a mortality weighting
        mort_w = 1

        #% CALCULATE IMPACT OF DENSITY DEPENDENT MORTAILTY----------#
        if ddm == true
            s_nhb = 0
            for i in nearby_ids(agent, model)
                if model[i].species_ID == species_ID
                    s_nhb += 1
                end
            end
            if s_nhb >= 4
                mort_w = (1 + ((s_nhb - 4) / 4) * 1)
            end
        end

        #% BASELINE MORTALITY---------------------------------------#  
        baseline_mortality = base_mortality[species_ID] + weather_adjustments

        #* Approximates the standard gap mortality model 
        #* (see Keane et al. 2001)
        if rand(Uniform(0.0, 1.0)) < (baseline_mortality .* mort_w)
            if gap_maker == 1
                expand[cell] = true
            end

            ## Record dying tree height as a cell list
            previous_height[cell] = a_height

            set_get_functions.kill_tree(
                id,
                model,
                grass,
                grass_invasion_prob,
                cell,
                )
            return

        #% MORTAILTY DUE TO AGE AND SUPPRESSION---------------------#
        #* trees older than 10 years have a chance of competition 
        #* based death
        elseif (
            age > 10.0 && 
            mean(previous_growth) < supp_tolerance[species_ID] &&
            rand(Uniform(0.0, 1.0)) < supp_mortality[species_ID]
            )
            if gap_maker == 1
                expand[cell] = true
            end
            
            ## Record dying tree height as a cell list
            previous_height[cell] = a_height

            set_get_functions.kill_tree(
                id,
                model,
                grass,
                grass_invasion_prob,
                cell,
                )
            return

        #% MORTALITY DUE TO PESTS----------------------------------#
        elseif (pests == true)
            pest_induced_mortality = rand(Normal(pest_mortality, pest_var))
            if rand(Uniform(0.0, 1.0)) < pest_induced_mortality
                if gap_maker == 1
                    expand[cell] = true
                end

                ## Record dying tree height as a cell list
                previous_height[cell] = a_height

                set_get_functions.kill_tree(
                    id,
                    model,
                    grass,
                    grass_invasion_prob,
                    cell,
                    )
                return
            end
        end
    end

    #//-------------------------------------------------------------------------------------------#
    #% EXPAND GAPS
    """
    # Expand gap around dying tree
    If the dying individual belongs to a gap-maker species then it creates a canopy gap potentially 
    extending beyond its own grid cell on its death. The tree is assumed to fall in a direction and
    potentially kill trees in that direction.
    ## Arguments:
    - `cell_ID::Int64`: Value matching the patch_id or position in model grid for target cell.
    - `model`: The AgentBasedModel object defining the current model. This object is usually 
    created using `Agents.ABM()`.
    - `grid::Matrix{Tuple{Int64, Int64}}`: Matrix containing all of the coordinates of cells in the
    current model.
    """
    function expand_gap(
        cell_ID::Int64,
        model,
        grid::Matrix{Tuple{Int64, Int64}},
        grass::Bool,
        grass_invasion_prob::Float64,
    )
        #* Calculate the number of cells over which the dying tree may fall
        h = trunc(Int64, (model.previous_height[cell_ID] / model.cell_grain) + 1)

        #% DETERMINE ALL POSSIBLY AFFECTED NEIGHBOURING TREES-------#
        nearby_trees = nearby_ids(grid[cell_ID], model::ABM{<:GridSpaceSingle}, h)
        ax = rand(1:2)
        fall_dir = model.pcor[cell_ID][1][ax]
        
        #% LOOP OVER ALL POSSIBLY AFFECTED TREES--------------------#
        for i in nearby_trees 
            #* If tree is in direction of fall record its height
            #* and species ID as a cell property and kill tree
            if model[i].pos[ax] == fall_dir
                cell = model[i].patch_here_ID
                model.previous_height[cell] = model[i].height
                
                set_get_functions.kill_tree(
                        model[i].id,
                        model,
                        grass,
                        grass_invasion_prob,
                        cell,
                        )
            end
        end
    end

    #//-------------------------------------------------------------------------------------------#
    #% CAPTURE GAPS
    """
    # Grow new tree in gap
    Potentially grow a new adult tree from the saplings bank for chosen empty cell. Each 
    sapling has a set annual probability of 25% of capturing a gap (i.e. transitioning to adult) if 
    there is no adult tree in the cell. Species ID is select via lottery competition weighted on 
    the basis of local light environment if there are multiple species. 
    ## Arguments:
    - `cell_ID::Int64`: Value matching the patch_id or position in model grid for target cell.
    - `model`: The AgentBasedModel object defining the current model. This object is usually 
    created using `Agents.ABM()`.
    - `seedlings::Vector{Vector{Int64}}`: Number of seedlings for each species for every cell in the
    model grid.
    - `saplings::Vector{Vector{Int64}}`: Number of saplings for each species for every cell in the
    model grid.
    - `nhb_light::Vector{Float64}`: Light environment for all cells in model grid. For more details
    regarding light environments see Dislich et *al* (2009)
    - `shade_tolerance::Vector{Float64}`: Species specific shade tolerances
    - `growth_forms::Vector{Int64}`: Growth type of agent, 1 = trees; 2 = tree ferns.
    - `b2_jabowas::Vector{Float64}`: List of species-specific allometric constant defined by Botkin 
    et *al* (1972) defined as `2 * (max_height - 1.37) / max_dbh` where 1.37 is the smallest adult 
    tree height.
    - `b3_jabowas::Vector{Float64}`: List of species-specific allometric constant defined by Botkin 
    et *al* (1972) defined as `(max_height - 1.37) / max_dbh ^ 2` where 1.37 is the smallest adult 
    tree height.
    - `last_change_tick::Vector{Int64}`: Value of previous tick in which the species ID in a cell
    changed for all cells in the model grid
    - `tick::Int64`: Current model timestep (i.e. tick).
    - `n_changes::Vector{Int64}`: Total number of species ID changes that have occured in a cell for
    all cells in the model grid
    """
    function capture_gap(
        cell_ID::Vector{Int64},
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
        new_agents_list
    )
        seedlings = seedlings[cell_ID]
        saplings = saplings[cell_ID]

        #% ASSESS WHETHER TO CAPTURE GAP----------------------------#
        #? 0.25 is probability of one sapling becoming an adult
        if sum(saplings[1]) >= 1 && rand(Uniform(0, 1)) < (1.0 - (0.25 ^ sum(saplings[1])))
            #* Create down weightings due to light environment and shade tolerance
            weights = floor.(Int, (abs.(nhb_light[cell_ID] .- 
                                        shade_tolerance)) * 100)

            regenbank_wgt = saplings[1] .* weights

            #% DEFINE NEW TREE/AGENT--------------------------------#
            new_species_id = distribution_functions.lottery(regenbank_wgt, true)
            new_species_id = new_species_id !== nothing ? new_species_id : rand(1:length(seedlings))
            #% For trees
            if growth_forms[new_species_id] == 1
                #* dbh (in m) of 0.01 m (1 cm) + noise (0, 0.01)
                dbh = 0.01 .+ rand(Uniform(0, 0.01))
                b2 = b2_jabowas[new_species_id]
                b3 = b3_jabowas[new_species_id]
                height = 1.37 .+ (b2 .* dbh) .- (b3 .* dbh .* dbh)
                age = 1.0
            #%For tree-ferns
            elseif growth_forms[new_species_id] == 2
                dbh = 0.01 .+ rand(Uniform(0, 0.01))
                height = 1.5 .+ rand(Uniform(0, 0.1))
                age = 1.0
            else
                error("The growth form for $new_species_id is undefined, please check species demography data")
            end

            last_change_tick[cell_ID] .= tick
            n_changes[cell_ID] .+= 1

            #% ASSIGN NEW TREE/AGENT TO GAP-------------------------#
            new_agent = [new_species_id,
                        cell_ID[1],
                        growth_forms[new_species_id],
                        height,
                        dbh,
                        age,]

            push!(new_agents_list, new_agent)

            #% EMPTY REGENERATION BANK------------------------------#
            seedlings = seedlings .- seedlings
            saplings = saplings .- saplings
        end
    end

    #//-------------------------------------------------------------------------------------------#
    #% REGENERATE SEEDLING BANK
    """
    # Sapling and seedling lifecycling
    Function undertakes mortality for the sapling and seedling banks in a cell and the growth of
    new saplings and seedlings. 
    ## Arguments:
    - `cell_ID::Int64`: Value matching the patch_id or position in model grid for target cell.
    - `seedlings::Vector{Vector{Int64}}`: Number of seedlings for each species for every cell in the
    model grid.
    - `saplings::Vector{Vector{Int64}}`: Number of saplings for each species for every cell in the
    model grid.
    - `seedling_mortality::Vector{Float64}`: Species specific baseline seedling mortality
    probabilities.
    - `sapling_mortality::Vector{Float64}`: Species specific baseline sapling mortality
    probabilities.
    - `seedling_transition::Vector{Float64}`: Species specific probability of a seedling 
    transitioning into a sapling in a tick.
    """
    function regenerate_patch_bank(
        cell_ID::Int64,
        seedlings::Vector{Vector{Int64}},
        saplings::Vector{Vector{Int64}},
        seedling_mortality::Vector{Float64},
        sapling_mortality::Vector{Float64},
        seedling_transition::Vector{Float64}
    )
        #% SEEDLING AND SAPLING MORTAILTY---------------------------#
        dead_saplings = rand.(Binomial.(abs.(saplings[cell_ID]), sapling_mortality))
        dead_seedlings = rand.(Binomial.(abs.(seedlings[cell_ID]), seedling_mortality))

        saplings[cell_ID] = saplings[cell_ID] - dead_saplings
        saplings[cell_ID] = replace(x-> x < 0 ? 0 : x, saplings[cell_ID])
        seedlings[cell_ID] = seedlings[cell_ID] - dead_seedlings
        seedlings[cell_ID] = replace(x-> x < 0 ? 0 : x, seedlings[cell_ID])

        #% ADD NEW SEEDLINGS AND SAPLINGS---------------------------#
        new_saplings = rand.(Binomial.(abs.(seedlings[cell_ID]), seedling_transition))

        saplings[cell_ID] = saplings[cell_ID] + new_saplings
        seedlings[cell_ID] = seedlings[cell_ID] + new_saplings
    end
end