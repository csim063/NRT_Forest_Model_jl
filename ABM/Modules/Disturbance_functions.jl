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
        model
    )
        if rand(Uniform(0,1)) < model.disturbance_freq
            #%Start disturbance
            grid = collect(positions(model))
            world_size = length(grid)
            disturbed_area = 0.0

            #Ask a random patch
            c = rand((1:world_size))
            model.disturbed[c] = true
            disturb_front = [c]
            max_area = rand(Exponential(model.max_disturb_size * world_size))

            #% Disturb spread
            #TODO  think this could be rewritten by just using patches with disturbed true??
            while length(disturb_front) > 0 && disturbed_area ≤ max_area
                new_disturb_front = []

                for i in disturb_front
                    for n in nearby_positions(grid[i], model::ABM{<:GridSpaceSingle}, 1)
                        #* 0.235 here is the approx Moore nhb percolation threshold for 
                        #* interesting bhv (see O'S & P 2013)
                        N = findall(x->x==[n], model.pcor)[1]
                        if model.disturbed[N] ≠ true && rand(Uniform(0,1)) ≤ 0.235
                            push!(new_disturb_front, N)
                            model.disturbed[N] = true
                        end
                    end
                end

                disturb_front = new_disturb_front
                disturbed_area += length(new_disturb_front)
            end

            #% ask disturbed patches a couple things
            for p in model.patch_ID[model.disturbed .== true]
                a_ID = id_in_position(grid[p], model::ABM{<:GridSpaceSingle})
                if a_ID ≠ 0
                    kill_agent!(a_ID, model)
                end

                n_species = length(model.seedlings[p])
                model.seedlings[p] = zeros(Int, n_species)
                model.saplings[p] = zeros(Int, n_species)
                model.disturbed[p] = false
            end

            #MOVE THIS INTO ABOVE FOR LOOP
            # n_species = length(model.seedlings[1])
            # model.seedlings[model.disturbed .== true][] .= zeros(Int, n_species)
            # model.saplings[model.disturbed .== true][] .= zeros(Int, n_species)
            # model.disturbed[model.disturbed .== true] .= false
        end
    end
end