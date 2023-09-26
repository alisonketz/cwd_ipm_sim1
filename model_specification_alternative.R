        for (k in 1:n_study_area) {
          for (t in 2:n_year) {

          #Female: set projection into population model matrix
          for (a in 2:n_agef - 1) {
            pop_sus[k, 1, a, t] <-  pop_sus[k, 1, a - 1, t - 1] *
                                        sn_sus[1, a - 1, t - 1] *
                                        (1 - psi[k, 1, a - 1, t - 1])
          }
          pop_sus[k, 1, n_agef, t] <- pop_sus[k, 1, n_agef - 1, t - 1] *
                                              sn_sus[1, n_agef - 1, t - 1] *
                                              (1 - psi[k, 1, n_agef - 1, t - 1]) +
                                              pop_sus[k, 1,  n_agef, t - 1] *
                                              sn_sus[1, n_agef, t - 1] *
                                              (1 - psi[k, 1, n_agef, t - 1])

          #Female: fawn class = total #females * unisex fawns per female/2
          pop_sus[k, 1, 1, t] <- (sum(pop_sus[k, 1, 1:n_agef, t - 1]) +
                                  sum(pop_inf[k, 1, 1:n_agef, t - 1])) *
                                  fec[t] * .5

          ##########
          ### Males
          ##########

          #Male: set projection into population model matrix
          for (a in 2:(n_agem - 1)) {
            pop_sus[k, 2, a, t] <- pop_sus[k, 2, a - 1, t - 1] *
                                  sn_sus[2, a - 1, t - 1] *
                                  (1 - psi[k, 2, a - 1, t - 1])
          }
          #Male: accumulating age class
          pop_sus_[k, 2, n_agem, t] <- pop_sus[k, 2, n_agem - 1, t - 1] *
                                      sn_sus[2, n_agem - 1, t - 1] *
                                      (1 - psi[k, 2, n_agem - 1, t - 1]) +
                                      pop_sus[k, 2,  n_agem, t - 1] *
                                      sn_sus[2, n_agem, t - 1] *
                                      (1 - psi[k, 2, n_agem, t - 1])

          # Male: fawn class = total #females * unisex fawns per female/2
          pop_sus[k, 2, 1, t] <- pop_sus[k, 1, 1, t]


          ###################################################
          ### Infected/Infectious
          ###################################################

          ###########
          ### Females
          ###########

          #Female: set projection into population model matrix
          for (a in 2:(n_agef - 1)) {
            pop_inf[k, 1, a, t] <- pop_inf[k, 1, a - 1, t - 1] *
                                        sn_inf[1, a - 1, t - 1] +
                                        pop_sus[k, 1, a - 1, t - 1] *
                                        sn_sus[1, a - 1, t - 1] *
                                        psi[k, 1, a - 1, t - 1]
          }
          ##Female: accumulating age = 9.5+ years
          pop_inf[k, 1, n_agef, t] <- pop_inf[k, 1, n_agef - 1, t - 1] *
                                      sn_inf[1, n_agef - 1, t - 1] +
                                      pop_inf[k, 1, n_agef, t - 1] *
                                      sn_inf[1, n_agef, t - 1] +
                                      pop_sus[k, 1, n_agef - 1, t - 1] *
                                      sn_sus[1, n_agef - 1, t - 1] *
                                      psi[k, 1, n_agef - 1, t - 1] +
                                      pop_sus[k, 1,  n_agef, t - 1] *
                                      sn_sus[1, n_agef, t - 1] *
                                      psi[k, 1, n_agef, t - 1] 

          ##Female: fawn class
          ##there are no infected fawns at birth
          pop_inf[k, 1, 1, t] <- 0

          ##########
          ### Males
          ##########

          ### Male: project forward anually
          for (a in 1:(n_agem - 1)) {
              pop_inf[k, 2, a, t] <- pop_inf[k, 2, a - 1, t - 1] *
                                            sn_inf[2, a - 1, t - 1] +
                                            pop_sus[k, 2, a - 1, t - 1] *
                                            sn_sus[2, a - 1, t - 1] *
                                            psi[k, 2, a - 1, t - 1]
          }

          ### Male: accumulating age class = 6.5+
          pop_inf[k, 2, n_agem, t] <- pop_inf[k, 2, n_agem - 1, t - 1] *
                                              sn_inf[2, n_agem - 1, t - 1] +
                                              pop_inf[k, 2, n_agem, t - 1] *
                                              sn_inf[2, n_agem, t - 1] +
                                              pop_sus[k, 2, n_agem - 1, t - 1] *
                                              sn_sus[2, n_agem - 1, t - 1] *
                                              psi[k, 2, n_agem - 1, t - 1] +
                                              pop_sus[k, 2,  n_agem, t - 1] *
                                              sn_sus[2, n_agem, t - 1] *
                                              psi[k, 2, n_agem, t - 1] 

          #Male: fawn class
          #there are no infected fawns at birth
          pop_inf[k, 2, 1, t] <- 0

          }
        }
