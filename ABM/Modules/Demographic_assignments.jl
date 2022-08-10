## TODO Documentation

module demog_metrics
    using Distributions

    function age_by_dbh(
        height,
        dbh,
        max_dbh,
        max_height,
        g_jabowa,
        b2_jabowa,
        b3_jabowa
    )
        ## Prepare counters for use in while loop
        # TODO Reassess this loop as it seems painfully complicated
        age = 0
        est_dbh = 0.01 + rand(Uniform(0, 0.01))[1]
        l_height = height

        ## Loop through until the desired dbh is reached
        while est_dbh ≤ dbh
            dbh_numerator = (est_dbh * g_jabowa * 
                                (1 - (est_dbh * l_height) /
                                    (max_dbh * max_height)))

            dbh_denomitor = (2.74 + 3 * b2_jabowa * est_dbh - 4 * b3_jabowa * 
                            est_dbh ^ 2)

            dbh_increment = dbh_numerator / dbh_denomitor

            l_height = 1.37 + (b2_jabowa * dbh) - (b3_jabowa * dbh * dbh)

            est_dbh += dbh_increment
            ## TODO may want to remove this hard coding
            age += rand(Normal(0.8, 0.1))
        end

        return age
    end

    function age_by_height(
        height,
        a_tf_height = 0.05289,
        b_tf_height = -0.05695
    )
        ## Prepare counters
        age = 0
        est_hgt = 0.0 + rand(Uniform(0, 0.01))[1]

        ## Loop through until the desired hgt is reached
        while est_hgt ≤ height
            hgt_increment = a_tf_height * exp(b_tf_height * est_hgt)

            est_hgt += hgt_increment
            age += rand(Normal(0.8, 0.1))
        end

        return age
    end
end