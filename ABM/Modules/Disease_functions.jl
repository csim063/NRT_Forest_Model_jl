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
                              transmission_age::Int64)
        #* Check whether agent gets infected by global chance
        #TODO add an increased chance if the agent is an edge patch
        if rand() < global_infection_prob
            agent.phytothera_infected::Bool = true
            
        #* If not globally infected check whether agent gets infected by local neighbours 
        #* if they are infected neighbours rand() < (1-(0.001 ^ 1))
        else
            for i in collect(nearby_ids(agent.pos, model::ABM, infectious_radius))
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
    function phytothera_impact(agent,)
        if agent.phytothera_infected == true
            agent.phytothera_infected_age += 1
        end

        #TODO Check infected above min symtomatic age and not already symptomatic then test to become symptomatic

        #TODO Check if symptomatic, if true then test to die
    end
end