# differential evolution optimazation approach
from numpy.random import rand
from numpy.random import choice
from numpy import asarray
from numpy import clip
from numpy import argmin
from numpy import min
from numpy import around
from Cycles.Cascade import cascade_cycle
import CoolProp as CP



# define objective function
def obj (refrigerant,parameters):
    #print(parameters)
    [COP,_] = cascade_cycle(refrigerant,parameters)
    return -COP



# define mutation operation
def mutation(x, F):
    return x[0] + F * (x[1] - x[2])


# define boundary check operation
def check_bounds(mutated, bounds):
    mutated_bound = [clip(mutated[i], bounds[i, 0], bounds[i, 1]) for i in range(len(bounds))]
    return mutated_bound


# define crossover operation
def crossover(mutated, target, dims, cr):
    # generate a uniform random value for every dimension
    p = rand(dims)
    # generate trial vector by binomial crossover
    trial = [mutated[i] if p[i] < cr else target[i] for i in range(dims)]
    return trial


def differential_evolution(refrigerant,bounds):
    # define population size
    pop_size = 20
    # define number of iterations
    iter = 50
    # define scale factor for mutation
    F = 0.5
    # define crossover rate for recombination
    cr = 0.7


    # initialise population of candidate solutions randomly within the specified bounds
    pop = bounds[:, 0] + (rand(pop_size, len(bounds)) * (bounds[:, 1] - bounds[:, 0])) # generate an aray of pop_size x len(X) with value from bound_min to bound_max
    print(pop)
    # evaluate initial population of candidate solutions
    cop_all = [obj(refrigerant,ind) for ind in pop] # cop_all is an array of pop_size x 1
    print('cop_all:',cop_all)
    # find the best performing vector of initial population
    best_vector = pop[argmin(cop_all)] #argmin Returns the indices of the minimum values along an axis.
    best_obj = min(cop_all)
    prev_obj = best_obj
    # run iterations of the algorithm
    for i in range(iter):
        # iterate over all candidate solutions
        for j in range(pop_size):
            # choose three candidates, a, b and c, that are not the current one
            candidates = [candidate for candidate in range(pop_size) if candidate != j]
            a, b, c = pop[choice(candidates, 3, replace=False)] # choose 3 random vector in pop, but not the current one
            # perform mutation
            mutated = mutation([a, b, c], F) # generate mutated vector from 3 randon vector
            #print('original mutated":', mutated)
            # check that lower and upper bounds are retained after mutation
            mutated = check_bounds(mutated, bounds)
            #print('mutated after check bound":', mutated)
            # perform crossover
            trial = crossover(mutated, pop[j], len(bounds), cr) # trial vector is either mutated vector or current vector in pop
            # compute objective function value for target vector
            print('trial',j,':',trial)
            obj_target = obj(refrigerant,pop[j])
            # compute objective function value for trial vector
            obj_trial = obj(refrigerant,trial)
            # perform selection
            if obj_trial < obj_target:
                # replace the target vector with the trial vector
                pop[j] = trial
                # store the new objective function value
                cop_all[j] = obj_trial
        print('iteration',i,':',cop_all)
        # find the best performing vector at each iteration
        best_obj = min(cop_all)
        # store the lowest objective function value
        if best_obj < prev_obj:
            best_vector = pop[argmin(cop_all)]
            prev_obj = best_obj
            # report progress at each iteration
            print('Iteration: ' ,i,'f(', around(best_vector, decimals=2),')', best_obj)
    return [best_vector, best_obj]

