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
end

module set_get_functions
    using Agents
    using StatsBase

    # TODO this function could be improved by vectorising
    function get_nhb_shade_height(
        cell_ID,
        model
    )
        grid = collect(positions(model))
        crit_heights = range(0, 36, step = 4)
        
        s_heights = Float64[]
        focal_height = model[id_in_position(grid[cell_ID], model::ABM{<:GridSpaceSingle})].height
        
        for n in range(1, length(model.nhb_set[cell_ID]))
            crit_height = crit_heights[n]
            for i in range(1, length(model.nhb_set[cell_ID][n]))
                pos_id = model.nhb_set[cell_ID][n][i]
                h = model[id_in_position(pos_id, model::ABM{<:GridSpaceSingle})].height
                if h > crit_height && h > focal_height
                        push!(s_heights, h)
                end
            end
        end
        
        return mean(s_heights)
    end

    function get_light_env(
        cell_ID,
        model,
        demography_df
    )
        l = model.nhb_shade_height[cell_ID] / maximum(demography_df.max_hgt)
        light = min(1.0, l)

        return light
    end

end
