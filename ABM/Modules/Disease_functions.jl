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
    """
    function phytothera_spread(agent,
                              model,
                              global_infection_prob::Float64,
                              phyto_local_prob::Float64,
                              infectious_radius::Int64)
        #* Check whether agent gets infected by global chance
        if rand() < global_infection_prob
            agent.phytothera_infected::Bool = true
        
        ## TODO: Maybe improve by creating a list of all infected agents and then only checking those
        ## TODO  for local infection
        #* If not globally infected check whether agent gets infected by local neighbours 
        #* if they are infected neighbours rand() < (1-(0.001 ^ 1))
        else
            for i in collect(nearby_ids(agent.pos, model::ABM, infectious_radius))
                if agent.phytothera_infected == true
                    break
                else model[i].phytothera_infected == true
                    if rand() < phyto_local_prob
                        agent.phytothera_infected::Bool = true
                    end
                end
            end
        end
    end
end