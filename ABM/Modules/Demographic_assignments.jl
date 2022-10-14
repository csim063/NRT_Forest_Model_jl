"""
Module contains functions to calculate the demographic parameters (currently just age) to be 
assigned to trees (agents) during initialisation either based on initial height for tree-ferns 
or dbh for trees.
"""
module demog_metrics
    using Distributions
    #//-------------------------------------------------------------------------------------------#
    #% CALCULATE AGENT AGE BASED ON DBH
    """
    # Age by DBH
    Use the diameter at breast height of a tree (agent) to calculate an initial age. This 
    calculation is based on the formula provided by Botkin et *al*. (1972). 
    ## References:
    - Botkin, D.B., Janak, J.F., Wallis, J.R., 1972. Some ecological consequences of acomputer model 
    of forest growth. J. Ecol. 60, 849-872.
    ## Arguments:
    - `height::Float64`: Height in meters of tree.
    - `dbh::Float64`: Diameter at breast height in meters of tree.
    - `max_dbh::Float64`: Species-specific maximum attainable diameter at breast height (m).
    - `max_height::Float64`: Species-specific maximum attainable height (m).
    - `g_jabowa::Float64`: Species-specific allometric constant defined by Botkin et *al* (1972).
    - `b2_jabowa::Float64`: Species-specific allometric constant defined by Botkin et *al* (1972)
        defined as `2 * (max_height - 1.37) / max_dbh` where 1.37 is the smallest adult tree height.
    - `b3_jabowa::Float64`: Species-specific allometric constant defined by Botkin et *al* (1972)
        defined as `(max_height - 1.37) / max_dbh ^ 2` where 1.37 is the smallest adult tree height.
    ## Return
    - Float64
    ## Examples
    ```julia-repl
    julia> age_by_dbh(agent.height, 
                    agent.dbh, 
                    model.max_dbhs[1], 
                    Float64(model.max_heights[1]), 
                    model.g_jabowas[1], 
                    model.b2_jabowas[1], 
                    model.b3_jabowas[1])
    74.90613803145392
    ```
    """
    function age_by_dbh(
        height::Float64,
        dbh::Float64,
        max_dbh::Float64,
        max_height::Float64,
        g_jabowa::Float64,
        b2_jabowa::Float64,
        b3_jabowa::Float64
    )
        #% DEFINE COUNTER VALUES USED IN WHILE LOOP-----------------#
        age = 0
        est_dbh = 0.01 + rand(Uniform(0, 0.01))
        l_height = height

        #% LOOP THROUGH UNTIL DESIRED DBH IS REACHED ---------------#
        #* Hard coded values are based on those used by 
        #* Morales and Perry (2017)
        while est_dbh ≤ dbh
            dbh_numerator = (est_dbh * g_jabowa * 
                                (1 - (est_dbh * l_height) /
                                    (max_dbh * max_height)))

            dbh_denomitor = (2.74 + 3 * b2_jabowa * est_dbh - 4 * b3_jabowa * 
                            est_dbh ^ 2)

            dbh_increment = dbh_numerator / dbh_denomitor

            l_height = 1.37 + (b2_jabowa * dbh) - (b3_jabowa * dbh * dbh)

            est_dbh += dbh_increment
            age += rand(Normal(0.8, 0.1))
        end

        return age
    end

    #//-------------------------------------------------------------------------------------------#
    #% CALCULATE AGENT AGE BASED ON HEIGHT
    """
    # Age by height
    Use the height a tree (agent) to calculate an initial age. This function is specifically 
    designed for use on agents with a growth form representing tree ferns. This calculation is based
    on the formula used by Brock et *al*. (2019). 
    ## References:
    - Brock, J. M. R., Morales, N. S., Burns, B. R., & Perry, G. L. W. (2019). The hare, tortoise 
    and crocodile revisited: Tree fern facilitation of conifer persistence and angiosperm growth in 
    simulated forests. The Journal of Ecology, 108(3), 969-981.
    - Botkin, D.B., Janak, J.F., Wallis, J.R., 1972. Some ecological consequences of acomputer model 
    of forest growth. J. Ecol. 60, 849-872.
    ## Arguments:
    - `height::Float64`: Height in meters of tree.
    - `a_tf_height::Float64`: Species-specific allometric constant defined by Brock et *al* (2019).
    *default = 0.05289*.
    - `b_tf_height::Float64`: Species-specific allometric constant defined by Brock et *al* (2019). 
    *default = -0.05695*.
    ## Return
    - Float64
    ## Examples
    ```julia-repl
    julia> age_by_height(agent.height)
    2.543424655317458
    """
    function age_by_height(
        height::Float64,
        a_tf_height::Float64 = 0.05289,
        b_tf_height::Float64 = -0.05695
    )
        #% DEFINE COUNTER VALUES USED IN WHILE LOOP-----------------#
        age = 0
        est_hgt = 0.0 + rand(Uniform(0, 0.01))[1]

        #% LOOP THROUGH UNTIL THE DESIRED HEIGHT IS REACHED---------#
        while est_hgt ≤ height
            hgt_increment = a_tf_height * exp(b_tf_height * est_hgt)

            est_hgt += hgt_increment
            age += rand(Normal(0.8, 0.1))
        end

        return age
    end
    
end

