# write function to take a greta model, simulate from priors and find valid free
# states as inits (those with all finite density and gradients)

# need to convert these outputs to inits
make_inits <- function(i, simulations, greta_arrays) {
    values <- lapply(simulations, greta.dynamics:::slice_first_dim, i)
    values <- mapply(enforce_dim, values, greta_arrays, SIMPLIFY = FALSE)
    do.call(initials, values)
}

# ensure the r value has the same dimensions as the corresponding greta array
enforce_dim <- function(r_value, greta_array) {
    array(r_value, dim = dim(greta_array))
}

# convert an inits list into a matrix of free states
get_free_states <- function(inits_list, variable_greta_arrays, model) {
    # attach the variable greta arrays here, to be found by
    # parse_initial_values(), which is hard-coded to look in a particular number
    # of parent frames above
    attach(variable_greta_arrays, warn.conflicts = FALSE)
    free_state_list <- lapply(inits_list,
                              greta:::parse_initial_values,
                              model$dag)
    do.call(rbind, free_state_list)
}

# sample a batch of n free state values from the model priors, but only return
# those that are valid
prior_sample_free_states_batch <- function(model, n) {

    # 1. simulate n times from priors for all parameters (as in calculate)
    variable_nodes <- model$dag$node_list[model$dag$node_types == "variable"]
    variable_greta_arrays <- lapply(variable_nodes, greta:::as.greta_array)
    sims <- do.call(calculate, c(variable_greta_arrays, list(nsim = n)))

    # 2. convert these back to free state values
    inits_list <- lapply(seq_len(n),
                         make_inits,
                         sims,
                         variable_greta_arrays)

    free_states <- get_free_states(inits_list,
                                   variable_greta_arrays,
                                   model)

    # 3. push the free state through to get density and gradients for all sims
    tfe <- model$dag$tf_environment
    # model$dag$define_joint_density()
    tfe$log_prob <- model$dag$generate_log_prob_function()
    tfe$joint_density_adj <- model$dag$on_graph(tfe$log_prob(tfe$free_state))
    tfe$grads <- tf$gradients(tfe$joint_density_adj, tfe$free_state)[[1]]
    tfe$density_grads <- tf$concat(list(tf$expand_dims(tfe$joint_density_adj, 1L),
                                        tfe$grads),
                                   axis = 1L)
    model$dag$send_parameters(free_states)
    density_grads <- tfe$sess$run(tfe$density_grads, feed_dict = tfe$feed_dict)

    # determine validity (finite density and grads) and return only the valid free states
    valid <- apply(is.finite(density_grads), 1, all)
    free_states[valid,  , drop = FALSE]

}

# Attempt to sample n valid free state values (those that result in a finite
# density and gradients) for the given model by sampling from the model priors.
# these can be use to define initial values for models that are difficult to
# sample from. Iteratively add samples to get at least n, but give up after
# trying max_tries.
prior_sample_free_states <- function(model, n, max_tries = n * 100, initial_tries = n) {

    # first attempt, hopefully they are all there
    free_states <- prior_sample_free_states_batch(model, initial_tries)
    tries <- initial_tries

    while (tries < max_tries & nrow(free_states) < n) {

        # calculate success rate and number still needed
        successes <- nrow(free_states)
        still_needed <- n - successes
        success_rate <- successes / tries

        # handle 0 success case
        finite_success_rate <- ifelse(success_rate == 0,
                                      1 / tries,
                                      success_rate)

        # work out how many to try to get there
        predicted_n_to_try <- still_needed / finite_success_rate
        # do some more, to increase the chance of getting there
        n_to_try <- round(1.1 * predicted_n_to_try)
        # don't try more than max_tries
        n_to_try <- pmin(n_to_try, max_tries - tries)

        # try them
        free_state_new <- prior_sample_free_states_batch(model, n_to_try)
        free_states <- rbind(free_states, free_state_new)
        tries <- tries + n_to_try

    }

    # see if we were successful
    successes <- nrow(free_states)

    # error informatively
    if (successes < n) {
        if (successes == 0) {
            msg <- cli::format_error(
                c("no valid initial values were found in {max_tries} samples",
                  "from the model priors")
            )
        } else {
            msg <- cli::format_error(
                c("only {successes} initial values were found in {max_tries} samples",
                  "from the model priors")
            )
        }

        stop(msg, call. = FALSE)

    }

    # otherwise just return the required number
    free_states[seq_len(n), ]

    # add a progress bar to this (progress =  successes/n)
    # maybe a second one simultaneously for tries / max_tries

}

# check that all the variable greta arrays are available so that inits can be
# specified for them
inits_are_deterministic <- function(model) {

    # check all the required greta arrays are visible
    visible_greta_arrays <- model$visible_greta_arrays
    visible_nodes <- lapply(visible_greta_arrays, greta:::get_node)
    visible_node_names <- vapply(visible_nodes,
                                 greta:::member,
                                 "unique_name",
                                 FUN.VALUE = character(1))

    variable_nodes <- model$dag$node_list[model$dag$node_types == "variable"]
    variable_node_names <- vapply(variable_nodes,
                                  greta:::member,
                                  "unique_name",
                                  FUN.VALUE = character(1))

    all(variable_node_names %in% visible_node_names)

}

# function to create a list of initial values object from a matrix of free states
initials_from_free_states <- function(model, free_states) {

    if(!inits_are_deterministic(model)) {
        msg <- cli::format_warning(
            c("not all variable greta arrays are visible in the current workspace",
              "some initial values for some variables cannot all be validated")
        )
        warning(msg, call. = FALSE)
    }

    n <- nrow(free_states)
    dag <- model$dag
    tfe <- dag$tf_environment
    visible_greta_arrays <- model$visible_greta_arrays

    # find visible greta arrays that are variables
    vga_nodes <- lapply(visible_greta_arrays, greta:::get_node)
    vga_node_names <- vapply(vga_nodes,
                             greta:::member,
                             "unique_name",
                             FUN.VALUE = character(1))
    all_variable_node_names <- names(dag$node_types[dag$node_types == "variable"])
    keep <- vga_node_names %in% all_variable_node_names
    visible_variable_greta_arrays <- visible_greta_arrays[keep]

    # hack the free states into calculate
    target_node_list <- lapply(visible_variable_greta_arrays, greta:::get_node)
    target_names_list <- lapply(target_node_list, dag$tf_name)
    target_tensor_list <- lapply(target_names_list, get, envir = tfe)
    assign("calculate_target_tensor_list", target_tensor_list, envir = tfe)
    dag$set_tf_data_list("batch_size", n)
    dag$build_feed_dict(list(free_state = free_states))

    sims <- dag$tf_sess_run("calculate_target_tensor_list", as_text = TRUE)

    # split the greta array values into inits with make_inits()
    inits_list <- lapply(seq_len(n),
                         make_inits,
                         sims,
                         visible_variable_greta_arrays)

    inits_list
}

# given a greta model, and a required number of chains, generate an initials
# object of the correct shape that ensures the initial values for all variable
# greta arrays visiblle in the workspace result in finite density and gradient
# values
generate_valid_inits <- function(model, chains, max_tries = chains * 100, initial_tries = chains) {
    free_states <- prior_sample_free_states(model,
                                            chains,
                                            max_tries = max_tries,
                                            initial_tries = initial_tries)
    initials_from_free_states(model,
                              free_states)
}
