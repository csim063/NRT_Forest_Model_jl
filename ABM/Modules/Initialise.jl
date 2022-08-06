"""
Module contains all of the functions required to setup the model at the start
of a run
"""
module Setup
    using Agents
    using Random
    using StatsBase
    #% Define agents
    #! !!!! BEST APPROACH IS PROBABLY TO HAVE ADULT TREES AS AGENTS AND 
    #! EVERYTHING ELSE INCLUDING SAPLINGS/SEEDLINGS AS PATCHE PROPERTIES
    #! Not that i have placed defining agents first as agents are required to
    #! actually initialise a model instance
    ## TODO: Currently just dummy agents to try get model working
    # TODO: probably will need to define different agents for different things
    #? Maybe a dict of species key value pairs could be useful
    @agent Tree GridAgent{2} begin 
        species_ID::Int
        height::Int
        dbh::Int
        age::Int
        treecolor::Symbol #Just used for
    end

    #% Define world
    function forest_model(;
        griddims = (50, 50),
        site_df = site_df,
        seed = 999,)
        
        #? Could we define space as many dimensional and just add the properties like height
        #? and species present to this
        space = GridSpaceSingle(griddims; periodic = false);
        rng = MersenneTwister(seed)

        model = ABM(Tree, space; 
            rng,
            scheduler = Schedulers.Randomly())

        ## Populate the world with adult tree agents
        grid = collect(positions(model))
        num_positions = prod(griddims)

        colour_list = [
            :firebrick1,
            :gold1,
            :turquoise3,
            :royalblue1,
            :deeppink2,
            :darkorange1,
            :forestgreen,
            :cadetblue2,
        ]

        #Make for loop that samples a proportion of space and allocates each species
        for p in 1:num_positions
            # Todo get correct heights etc for each tree
            #? Could we use dictionary keys to get name value pairs and make it clearer what we are doing
            # Column 1 is species column 2 is initial abundance
            specID = wsample(site_df[ : , 1], site_df[ : , 2])

            adult_tree = Tree(
                p,
                grid[p],
                specID,
                0,
                0,
                0,
                colour_list[specID]
            )
            add_agent_single!(adult_tree, model)
        end
        
        return model
    end

    #% Define patches

end