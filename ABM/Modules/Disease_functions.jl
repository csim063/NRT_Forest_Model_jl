"""
Module contains all functions related to the spread and impact of the different disease types
implemented in the model.
"""

module disease_functions
    using Agents

    #//-------------------------------------------------------------------------------------------#
    #% PHYTOTHERA SPREAD
    """
    Function implements the spread of phytothera disease through soils in the model. This disease
    TODO: Add more detail here use some of the details you have in obsidian
    TODO: Refactor to avoid need to access globals agent and model
    """
    function phytothera_spread(agent,
                              model,
                              global_infection_prob::Float64,
                              phyto_local_prob::Float64,
                              infectious_radius::Int64,
                              transmission_age::Int64,
                              pos::Tuple{Int64, Int64},)
        #* Check whether agent gets infected by global chance
        #TODO add an increased chance if the agent is an edge patch
        if rand() < global_infection_prob
            agent.phytothera_infected::Bool = true
            
        #* If not globally infected check whether agent gets infected by local neighbours 
        #* if they are infected neighbours rand() < (1-(0.001 ^ 1))
        else
            for i in collect(nearby_ids(pos, model::ABM, infectious_radius))
                if agent.phytothera_infected == true
                    break
                else model[i].phytothera_infected_age ≥ transmission_age && rand() < phyto_local_prob
                        agent.phytothera_infected::Bool = true
                end
            end
        end
    end

    #//-------------------------------------------------------------------------------------------#
    #% PHYTOTHERA IMPACT
    """
    TODO: Document this function
    TODO: Refactor to avoid need to access globals agent and model
    """
    function phytothera_impact(agent,
                               model,
                               phytothera_infected::Bool,
                               min_symptomic_age::Int64,
                               symptom_prob::Float64,
                               mortality_prob::Float64,
                               gap_maker::Int64,
                               species_ID::Int64,
                               cell::Int64,
                               expand::BitVector,
                               previous_species::Vector{Float64},
                               previous_height::Vector{Float64},
                               a_height::Float64,
                               id::Int64,
                               )

        if phytothera_infected == true
            agent.phytothera_infected_age += 1
        end

        #* Check whether tree has become symptomatic
        if agent.phytothera_infected_age ≥ min_symptomic_age && 
            agent.phytothera_symptomatic == false
            if rand() < symptom_prob
                agent.phytothera_symptomatic::Bool = true
            end
        end

        #* If symptomatic chance for tree to die
        if agent.phytothera_symptomatic == true 
            if rand() < mortality_prob
                if gap_maker[species_ID] == 1
                    expand[cell] = true
                end
    
                ## Record dying tree species and height as a cell list
                previous_species[cell] = species_ID
                previous_height[cell] = a_height
    
                kill_agent!(id, model)
            end    
        end
    end
end