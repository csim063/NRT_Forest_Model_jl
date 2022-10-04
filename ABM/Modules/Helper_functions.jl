## TODO Write a description

module distribution_functions
    using Distributions
    ## TODO Document
    function generate_LogNormal(m,std)
        γ = 1+std^2/m^2
        μ = log(m/sqrt(γ))
        σ = sqrt(log(γ))

        return LogNormal(μ,σ)
    end

    """
    Function to allocate species based on saplings
    Note it combines a number of Netlogo functions into one
    """
    function lottery(
        input_list,
        cumulative = true
        )

        summed_list = sum(input_list)[1]

        if summed_list > 0
            scale_bank = input_list ./ summed_list
            if cumulative == true
                scale_bank = cumsum(scale_bank)
            end

            #! I do not think this approach works if not set to cumulative above
            #! Also that findfirst syntax with [] and . is nuts
            r = rand(Uniform(0,1))
            r_indx = findfirst(x->x[1] .> r[1], scale_bank)
            
            return(r_indx[1])
        end
    end
end

"""
TODO
"""
module set_get_functions
    using Agents
    using StatsBase

    # TODO this function could be improved by vectorising
    function get_nhb_shade_height(
        cell_ID,
        model,
        grid,
        crit_heights,
        shell_width
    )
        if id_in_position(grid[cell_ID], model::ABM{<:GridSpaceSingle}) == 0
            final_s_hgts = 0.0
        else
            s_heights = Float64[]
            agent = model[id_in_position(grid[cell_ID], model::ABM{<:GridSpaceSingle})]
            focal_height = agent.height

            for n in 1:shell_width
                if n ≠ 1
                    for idx in setdiff(nearby_ids(agent.pos, model, n), nearby_ids(agent.pos, model, n-1))
                        if model[idx].height > crit_heights[n] && model[idx].height > focal_height
                            push!(s_heights, model[idx].height)
                        end
                    end
                else
                    for idx in nearby_ids(agent.pos, model, 1)
                        if model[idx].height > crit_heights[n] && model[idx].height > focal_height
                            push!(s_heights, model[idx].height)
                        end
                    end
                end
            end
            final_s_hgts = length(s_heights) > 0 ? mean(s_heights) : 0.0
        end
        return final_s_hgts
    end

    function get_light_env(
        nhb_shade_height::Float64,
        max_heights::Vector{Int64}
    )
        l = nhb_shade_height ./ maximum(max_heights)
        light = min(1.0, l)

        return light
    end

    """
    Simple function to update count of total trees of each species in the model
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

    """
    Function to select random neighbor and increase seedlings there
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
end
