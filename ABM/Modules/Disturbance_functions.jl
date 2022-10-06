## TODO Write a description

module disturbance_functions
    using Agents
    using StatsBase
    using Distributions
    using FillArrays

    """ 
    Function to implement landscape level disturbances
    """
    #TODO remember that currently you are calling this as a model step
    #TODO run only once per tick
    function lsp_disturbance(
        model,
        grid::Matrix{Tuple{Int64, Int64}},
        disturbance_freq::Float64,
        disturbed::BitVector,
        max_disturb_size::Float64,
        #nhb_set::Vector{Vector{Tuple{Int64, Int64}}},
        nhb_set_ids::Vector{Vector{Int64}},
        close_nhbs_count::Vector{Int64},
        pcor::Vector{Vector{Tuple{Int64, Int64}}},
        p_ID::Vector{Int64},
        n_species::Int64,
        seedlings::Vector{Vector{Int64}},
        saplings::Vector{Vector{Int64}}
    )
        if rand(Uniform(0,1)) < disturbance_freq
            #%Start disturbance
            world_size = length(grid)
            disturbed_area = 0.0

            #Ask a random patch
            c = rand((1:world_size))
            disturbed[c] = true
            disturb_front = [c]
            max_area = rand(Exponential(max_disturb_size * world_size))

            #% Disturb spread
            #TODO  think this could be rewritten by just using patches with disturbed true??
            while length(disturb_front) > 0 && disturbed_area ≤ max_area
                new_disturb_front = []

                for i in disturb_front
                    #for n in nearby_positions(grid[i], model::ABM{<:GridSpaceSingle}, 1)
                    for N in nhb_set_ids[i][1:close_nhbs_count[i]]
                        #* 0.235 here is the approx Moore nhb percolation threshold for 
                        #* interesting bhv (see O'S & P 2013)
                        #N = findfirst(isequal([n]), pcor)
                        if disturbed[N] ≠ true && rand(Uniform(0,1)) ≤ 0.235
                            push!(new_disturb_front, N)
                            disturbed[N] = true
                        end
                    end
                end

                disturb_front = new_disturb_front
                disturbed_area += length(new_disturb_front)
            end

            #% ask disturbed patches a couple things
            for p in p_ID[disturbed .== true]
                a_ID = id_in_position(grid[p], model::ABM{<:GridSpaceSingle})
                if a_ID ≠ 0
                    kill_agent!(a_ID, model)
                end

                seedlings[p] = zeros(Int, n_species)
                saplings[p] = zeros(Int, n_species)
                disturbed[p] = false
            end
        end
    end
end