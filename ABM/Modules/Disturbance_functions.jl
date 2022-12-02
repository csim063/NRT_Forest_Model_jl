"""
Module contains all functions which are used to create disturbances in the model.
"""
module disturbance_functions
    using Agents
    using StatsBase
    using Distributions
    using FillArrays

    #//-------------------------------------------------------------------------------------------#
    #% LANDSCAPE LEVEL DISTURBANCE
    """ 
    # Landscape disturbances
    Function implements landscape disturbances which spread via percolation. This could be 
    considered a disturbance such as fire.
    ## Arguments:
    - `model`: The AgentBasedModel object defining the current model. This object is usually 
    created using `Agents.ABM()`.
    - `grid::Matrix{Tuple{Int64, Int64}}`: Matrix containing all of the coordinates of cells in the
    current model.
    - `disturbance_freq::Float64`: Probability of a disturbance happening in any single tick
    - `disturbed::BitVector`: Vector for every cell in the model indicating whether the cell is 
    considered  currently disturbed (1) or undisturbed (0).
    - `max_disturb_size::Float64`: Maximum proportion (0-1) of model cells that may be disturbed at
    a single time.
    - `nhb_set_ids::Vector{Vector{Int64}}`: Patch_IDs of all cells considered neighbours (i.e.
    considered to effect) for all cells in the model grid.
    - `close_nhbs_count::Vector{Int64}`: Number of neighbouring cells in distance one of each cell
     in the model grid.
    - `p_ID::Vector{Int64}`: Patch_IDs for all patches in model grid.
    - `n_species::Int64`:  Number of total species.
    - `seedlings::Vector{Vector{Int64}}`: Number of seedlings for each species for every cell in the
    model grid.
    - `saplings::Vector{Vector{Int64}}`: Number of saplings for each species for every cell in the
    model grid.
    ```
    """
    function lsp_disturbance(
        model,
        grid::Matrix{Tuple{Int64, Int64}},
        disturbance_freq::Float64,
        disturbed::BitVector,
        max_disturb_size::Float64,
        nhb_set_ids::Vector{Vector{Int64}},
        close_nhbs_count::Vector{Int64},
        p_ID::Vector{Int64},
        n_species::Int64,
        seedlings::Vector{Vector{Int64}},
        saplings::Vector{Vector{Int64}},
        grass::Bool,
        grass_invasion_prob::Float64,
    )
        ## Check if a disturbance has occurred this tick
        if rand(Uniform(0,1)) < disturbance_freq
            #% START DISTURBANCE------------------------------------#
            world_size = length(grid)
            disturbed_area = 0.0

            #* Select a random patch to be the start of the disturbance front
            c = rand((1:world_size))
            disturbed[c] = true
            disturb_front = [c]
            max_area = rand(Exponential(max_disturb_size * world_size))

            #% SPREAD DISTURBANCE-----------------------------------#
            while length(disturb_front) > 0 && disturbed_area ≤ max_area
                new_disturb_front = []

                for i in disturb_front
                    for N in nhb_set_ids[i][1:close_nhbs_count[i]]
                        #? 0.235 here is the approx Moore nhb percolation threshold for 
                        #?  (see O'Sullivan & Perry 2013)
                        if disturbed[N] ≠ true && rand(Uniform(0,1)) ≤ 0.235
                            push!(new_disturb_front, N)
                            disturbed[N] = true
                        end
                    end
                end

                disturb_front = new_disturb_front
                disturbed_area += length(new_disturb_front)
            end

            #% SET VALUES FOR DISTURBED PATCHES---------------------#
            for p in p_ID[disturbed .== true]
                a_ID = id_in_position(grid[p], model::ABM{<:GridSpaceSingle})
                if a_ID ≠ 0
                    set_get_functions.kill_tree(
                        a_ID,
                        model,
                        grass,
                        grass_invasion_prob,
                        cell,
                    )
                end

                #* Reset seedlings and saplings to 0
                seedlings[p] = zeros(Int, n_species)
                saplings[p] = zeros(Int, n_species)
                disturbed[p] = false
            end
        end
    end
end