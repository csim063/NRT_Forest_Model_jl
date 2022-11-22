"""
Script contains two modules both of which contain functions which are largely miscellaneous and do 
not easily fit into other modules, however are helpful for model scripting. 

The first module provides functions which define distributions not easily created in Julia as well
as functions to streamline making use of randomness.

The second module contains functions to create agent sets, assign values to agent in agent sets, and
report values from agent sets.
"""

"""
Module provides functions which define distributions not easily created in Julia as wellas functions 
to streamline making use of randomness.
"""
module distribution_functions
    using Distributions

    #//-------------------------------------------------------------------------------------------#
    #% GENERATE A VALUE FROM A LOG NORMAL DISTRIBUTION
    """
    # Generate Log Normal Distribution
    Use the diameter at breast height of a tree (agent) to calculate an initial age. This 
    calculation is based on the formula provided by Botkin et *al*. (1972). 
    ## Arguments:
    - `m::Float64`: Mean of desired distribution.
    - `std::Float64`: Standard deviation of desired distribution
    ## Return
    - LogNormal{Float64}
    ## Examples
    ```julia-repl
    julia> generate_LogNormal(1.0, 0.5)
    LogNormal{Float64}(μ=-0.11157177565710491, σ=0.47238072707743883)
    ## To access a random value from distribution
    julia> rand(generate_LogNormal(1.0, 0.5))
    0.8114035723877778
    ```
    """
    function generate_LogNormal(m::Float64,
                                std::Float64)
        γ = 1.0+std^2.0/m^2.0
        μ = log(m/sqrt(γ))
        σ = sqrt(log(γ))

        return LogNormal(μ,σ)
    end

    #//-------------------------------------------------------------------------------------------#
    #% ASSIGN SPECIES IDS
    """
    # Lottery assign
    Function to allocate random species ID weighted by the number of saplings in a patch.
    Note it combines a number of Netlogo functions into one
    ## Arguments:
    - `input_list::Vector{Int64}`: Weighted list of number of saplings of each species.
    - `cumulative::Bool`: Whether to use cumulative sum list during assignment. *default = true*
    ## Return
    - Int64, describing species ID based on positions supplied in input_list
    ## Examples
    ```julia-repl
    julia> lottery([1, 1, 1, 1, 0, 0, 0, 1])
    8
    ```
    """
    function lottery(
        input_list::Vector{Int64},
        cumulative::Bool = true
        )

        summed_list = sum(input_list)[1]

        if summed_list > 0
            scale_bank = input_list ./ summed_list
            if cumulative == true
                scale_bank = cumsum(scale_bank)
            end

            r = rand(Uniform(0,1))
            r_indx = findfirst(x->x[1] .> r[1], scale_bank)
            
            return(r_indx[1])
        end
    end
end

#//------------------------------------------------------------------------------------------------#
#%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#//------------------------------------------------------------------------------------------------#

"""
Module contains functions to create agent sets, assign values to agent in agent sets, and report 
values from agent sets.
"""
module set_get_functions
    using Agents
    using StatsBase

    #//-------------------------------------------------------------------------------------------#
    #% GET HEIGHTS OF NEIGHBOURING TREES
    """
    # Get neighbourhood shade height
    Function obtains the mean height of neighbouring trees with height greater than a critical 
    value. The slected trees are distance weighted such that cells in distance i have to be of a
    height greater than item i in the critical heights list.
    ## Arguments:
    - `cell_ID::Int64`: Value matching the patch_id or position in model grid for target cell. Note
    that target cell is the central cell whose neighbourhood is going to be calculated.
    - `model`: The AgentBasedModel object defining the current model. This object is usually 
    created using `Agents.ABM()`.
    - `grid::Matrix{Tuple{Int64, Int64}}`: Matrix containing all of the coordinates of cells in the
    current model.
    - `crit_heights`: Iterator matching the critical height values weighting each distance from the
    target cell.
    - `shell_width::Int64`: Maximum distance value over which to obtain neighbour heights. I.E. if
    `shell_width = 2` the neighbour heights for all trees within a distance of 2 cells and greater
    than critical height for that distance will be obtained.
    ## Return
    - Float64
    ## Examples
    ```julia-repl
    julia> get_nhb_shade_height(1, 
                                model,
                                collect(positions(model)),
                                range(0, 32, step = 4),
                                8)
    7.8730749204029715
    ```
    """
    function get_nhb_shade_height(
        cell_ID::Int64,
        model,
        grid::Matrix{Tuple{Int64, Int64}},
        crit_heights,
        shell_width::Int64
    )
        #* Only calculate neighbourhood heights if there is a tree in the current cell
        if id_in_position(grid[cell_ID], model::ABM{<:GridSpaceSingle}) == 0
            final_s_hgts = 0.0
        else
            #% DEFINE COUNTER VALUES USED IN FOR LOOP---------------#
            s_heights = Float64[]
            agent = model[id_in_position(grid[cell_ID], model::ABM{<:GridSpaceSingle})]
            focal_height = agent.height
            focal_pos = agent.pos

            #% ITERATE OVER EACH CHOSEN DISTANCE FROM TARGET CELL---#
            for n in getindex(shell_width)
                ## If distance greater than one find only the targets at exactly that distance
                ## TODO try ∉ approach to remove setdiff
                if n ≠ 1
                    for idx in setdiff(nearby_ids(focal_pos, model, n), 
                                    nearby_ids(focal_pos, model, n-1))
                        idx_height = model[idx].height
                        if idx_height > crit_heights[n] && idx_height > focal_height
                            push!(s_heights, idx_height)
                        end
                    end
                ## If distance one just look for agents at distance one (more efficient)
                else
                    for idx in nearby_ids(focal_pos, model, 1)
                        idx_height = model[idx].height
                        if idx_height > crit_heights[n] && idx_height > focal_height
                            push!(s_heights, idx_height)
                        end
                    end
                end
            end

            final_s_hgts = length(s_heights) > 0 ? mean(s_heights) : 0.0
        end

        return final_s_hgts
    end

    #//-------------------------------------------------------------------------------------------#
    #% CALCULATE LIGHT ENVIRONMENT
    """
    # Get light/shade environment
    Based on the mean neighbourhood shading height (see `get_nhb_shade_height()`) calculate the 
    light environment. This light environment impacts the growth and gap capture of the tree 
    species. This approach is adopted from Dislich et *al* (2009). Note that 0 is high an 1 is low 
    in this approach.
    ## References:
    - Dislich, C., Günter, S., Homeier, J., Schröder, B., & Huth, A. (2009). Simulating Forest 
    Dynamics of a Tropical Montane Forest in South Ecuador. Erdkunde, 63(4), 347-364.
    ## Arguments:
    - `nhb_shade_height::Float64`: Mean weighted height of neighbouring trees as calculated by
    `get_nhb_shade_height()`.
    - `max_heights::Vector{Int64}`: Species-specific maximum attainable height (m). 
    ## Return
    - Float64, scaled between 0 and 1
    ## Examples
    ```julia-repl
    julia> get_light_env(model.nhb_shade_height[1], 
                        model.max_heights)
    0.21016138020700623
    ```
    """
    function get_light_env(
        nhb_shade_height::Float64,
        max_heights::Vector{Int64}
    )
        l = nhb_shade_height ./ maximum(max_heights) # Max max achievable height accross all species
        light = min(1.0, l)

        return light
    end

    #//-------------------------------------------------------------------------------------------#
    #% CALCULATE TOTAL TREES
    """
    # Update species abundances
    Simple function to update count of total trees of each species in the model
    ## Arguments:
    - `model`: The AgentBasedModel object defining the current model. This object is usually 
    created using `Agents.ABM()`.
    - `n_species::Int64`: Number of total species.
    ## Examples
    ```julia-repl
    julia> update_abundances(model, 8)
    ```
    """
    function update_abundances(
        model,
        n_species::Int64
    )
        model.abundances = zeros(Int64, n_species)

        for i in allids(model)
            model.abundances[model[i].species_ID] += Int64(1)
        end
    end

    #//-------------------------------------------------------------------------------------------#
    #% DISPERSE SEED TO RANDOM NEIGHBOUR
    """
    # Assign seedlings
    Function assign *N* seedlings randomly to neighbouring cells. Note that 1 seedling is assigned 
    at a time so each cell may receive multiple seedlings and/or multiple cells may receive 
    seedlings
    ## Arguments:
    - `n_seeds::Int64`: Number of seeds to assign
    - `dispersal_nhb_id::Vector{Int64}`: Patch_IDs for all neighbouring patches over which seeds 
    may be dispersed
    - `shad_h::Vector{Float64}`: Mean shade heights as calculated by `get_nhb_shade_height()` for 
    all cells in model grid
    - `r_hgt::Int64`: Regeneration height in meters for selected species
    - `seedlings::Vector{Vector{Int64}}`: Number of seedlings for each species for every cell in the
    model grid.
    - `spec_ID::Int64`: Selected species ID.
    ## Examples
    ```julia-repl
    julia> assign_seedling(23,
                        model.nhb_set_ids[agent.patch_here_ID]1:8],
                        model.nhb_shade_height,
                        35,
                        model.seedlings,
                        8)
    ```
    """
    function assign_seedling(
        n_seeds::Int64,
        dispersal_nhb_id::Vector{Int64},
        shad_h::Vector{Float64},
        r_hgt::Int64,
        seedlings::Vector{Vector{Int64}},
        spec_ID::Int64
    )
        for _ in 1:n_seeds
            rand_cell_ID = rand(dispersal_nhb_id)
            
            if shad_h[rand_cell_ID] ≤ r_hgt
                seedlings[rand_cell_ID][spec_ID] += Int64(1)
            end
        end
    end

    #//-------------------------------------------------------------------------------------------#
    #% GET NEIGHBOURS AT SPECIFIC DISTANCE
    #TODO CORRECT DOCUMENTAION
    """
    # Get neighbours at specific distance
    Function to get the patch_IDs of all neighbours at a specific distance from the focal cell.
    ## Arguments:
    - `pos::Tuple{Int64, Int64}`: Position of focal cell
    - `model`: The AgentBasedModel object defining the current model. This object is usually
        created using `Agents.ABM()`.
    - `D::Int64`: Distance from focal cell to get neighbours
    ## Return
    - Vector{Tuple{Int64, Int64}}, patch positions of all neighbours at distance `D` from focal cell
    ## Examples
    ```julia-repl
    julia> positions_at((1,1), model, 1)
    [(1, 2), (2, 1), (2, 2)]
    ```
    """ 
    function positions_at(
        pos::Tuple{Int64, Int64},
        model::ABM{<:GridSpaceSingle},
        D::Int64
    )
        #TODO: Try use nearby_positions() source code altered to only return positions at distance D

        a = collect(nearby_positions(pos::Tuple{Int64, Int64}, 
                                    model::ABM{<:GridSpaceSingle}, 
                                    (D::Int64-1)));

        D_nhbs = collect(nearby_positions(pos::Tuple{Int64, Int64}, 
                                        model::ABM{<:GridSpaceSingle}, 
                                        D::Int64));
        setdiff!(D_nhbs::Vector{Tuple{Int64, Int64}}, a::Vector{Tuple{Int64, Int64}})

        return D_nhbs::Vector{Tuple{Int64, Int64}}
    end
end

