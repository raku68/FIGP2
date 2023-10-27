import time
import random
import statistics
from matplotlib.pyplot import flag
import copy
import itertools

import numpy as np
import pandas as pd
from deap import tools
from deap.algorithms import varAnd, varOr

from .log_manager import table_log, txt_log

def unique_varAnd(population, toolbox, cxpb, mutpb,
                  check_func=None, max_trial=20):
    
    offspring = [toolbox.clone(ind) for ind in population]

    # Apply crossover and mutation on the offspring
    for i in range(1, len(offspring), 2):
        if random.random() < cxpb:
            olds = [str(offspring[i - 1]), str(offspring[i])]
            
            parent_1 = toolbox.clone(offspring[i - 1])
            parent_2 =  toolbox.clone(offspring[i])
            
            Decision_offsp_1, Decision_offsp_2 = False, False
            
            for _ in range(max_trial):
                candidate4offspring_1, candidate4offspring_2 = toolbox.mate(toolbox.clone(parent_1), toolbox.clone(parent_2))
                
                if (Decision_offsp_1 == False & check_func.is_unobserved(candidate4offspring_1, add_pool=False)):
                    offspring[i - 1] = candidate4offspring_1
                    Decision_offsp_1 = True
                    
                if (Decision_offsp_2 == False & check_func.is_unobserved(candidate4offspring_2, add_pool=False)):
                    offspring[i] = candidate4offspring_2
                    Decision_offsp_2 = True
                
                if (Decision_offsp_1 & Decision_offsp_2):
                    break
                
            if Decision_offsp_1 == False:
                offspring[i - 1] = candidate4offspring_1
            if Decision_offsp_2 == False:
                offspring[i] = candidate4offspring_2

            del offspring[i - 1].fitness.values, offspring[i].fitness.values

    for i in range(len(offspring)):
        if random.random() < mutpb:
            parent_1 = toolbox.clone(offspring[i])
            for _ in range(max_trial):
                candidate4offspring_1, = toolbox.mutate(toolbox.clone(parent_1))
                if check_func.is_unobserved(candidate4offspring_1, add_pool=False):
                    break
            offspring[i] = candidate4offspring_1
            del offspring[i].fitness.values

    return offspring


def FGP_NLS_algorithm(population, toolbox, cxpb, mutpb, ngen,# stats=None,
             halloffame=None, num_elite_select=1,
             var_max_trial = 20, check_func=None,
             text_log=None, save_dir='./results', 
             func_name=['add', 'sub', 'mul', 'div', 'ln', 'exp', 'sqrt', 'square', 'cube'],
             n_expr2save=3
             ):
    
    xname = toolbox.compile.keywords['pset'].arguments
    pop_analysis_pool, pop_analysis_pool_select = list(), list()
    tab_log = table_log()

    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in population if not ind.fitness.valid]
    
    text_log.print(['=============== gen 0 ==============='])
    text_log.print([f'ind {_en} : {_ind}' for _en, (_ind, _) in enumerate(zip(invalid_ind, range(n_expr2save)))])
    
    fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit

    if halloffame is not None:
        halloffame.update(population)
    
    fitnesses = [ind.fitness.values[0] for ind in invalid_ind]
    tab_log.save_log( gen=0,
                      score_list      = [min(fitnesses), 
                                         statistics.median(fitnesses),
                                         max(fitnesses),
                                         np.nanstd(np.where(np.isinf(fitnesses), np.nan, fitnesses)),
                                         len(set([str(ind) for ind in invalid_ind]))/len(fitnesses),
                                         len(invalid_ind)],
                       
                      score_name_list = ['score-min', 'score-med', 'score-max', 'score-std', 'unique_rate', 'nevals'])
    pop_analysis_pool.append(pop_analysis(population, func_name=func_name, x_name=xname))

    # Begin the generational process
    for gen in range(1, ngen + 1):
        _start_gen = time.time()
        
        # Select the next generation individuals
        offspring = toolbox.select(population, len(population))
        
        pop_analysis_pool_select.append(pop_analysis(offspring, func_name=func_name, x_name=xname))
        _now_gen_log = [[f'=============== gen {gen} ===============']]
        _now_gen_log.extend([[f'ind {_en}', f'After selection : {str(_ind)}\t\t{_ind.fitness.values}'] for _en, (_ind, _) in enumerate(zip(offspring, range(n_expr2save)))]) 

        # Vary the pool of individuals
        if var_max_trial != 1:
            offspring = unique_varAnd(offspring, toolbox, cxpb, mutpb,
                                      check_func=check_func, max_trial=20)
        else:
            offspring = varAnd(offspring, toolbox, cxpb, mutpb)
        
        [_now_gen_log[_en+1].extend([f'After evolution : {str(_ind)}\t\t{_ind.fitness.values}']) for _en, (_ind, _) in enumerate(zip(offspring, range(n_expr2save)))]

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        [_now_gen_log[_en+1].extend([f'After const opt : {str(_ind)}\t\t{_ind.fitness.values}\nexpr state\n={_ind.state}']) for _en, (_ind, _) in enumerate(zip(offspring, range(n_expr2save)))]

        _now_gen_log = list(itertools.chain.from_iterable(_now_gen_log))
        text_log.print(_now_gen_log)
        
        # Update the hall of fame with the generated individuals
        if halloffame is not None:
            halloffame.update(offspring)
        
        # Replace the current population by the offspring
        population[:] = offspring
        
        for ind in offspring:
            check_func.is_unobserved(ind, add_pool=True)
        
        if num_elite_select > 0:
            for i in range(num_elite_select):
                population[-(i+1)] = halloffame[i]
        
        fitnesses = [ind.fitness.values[0] for ind in population]
        log = tab_log.save_log( gen=gen,
                                score_list = [min(fitnesses),
                                              statistics.median(fitnesses),
                                              max(fitnesses),
                                              np.nanstd(np.where(np.isinf(fitnesses), np.nan, fitnesses)),
                                              len(set([str(ind) for ind in invalid_ind]))/len(fitnesses),
                                              len(invalid_ind)], 
                                score_name_list = ['score-min', 'score-med', 'score-max', 'score-std', 'unique_rate', 'nevals'])
        pop_analysis_pool.append(pop_analysis(population, func_name=func_name, x_name=xname))
        text_log.print([
            f'best ind       : {str(halloffame[0])}\nbest score     : {halloffame[0].fitness.values[0]}\nExecution time : {time.time() - _start_gen:.4f} s', 
            f'=============== gen {gen} ===============\n']
            )
    
    pop_analysis_pool = pd.DataFrame(pop_analysis_pool)
    pop_analysis_pool.to_csv(f'{save_dir}/001_GP_node_analysis.tsv', sep='\t')
    pop_analysis_pool_select = pd.DataFrame(pop_analysis_pool_select)
    pop_analysis_pool_select.to_csv(f'{save_dir}/001_GP_node_analysis_select.tsv', sep='\t')

    return population, log

def isfloat(_str):
    try:
        float(_str)
    except:
        return False
    return True

def pop_analysis(pop, func_name, x_name):
    def count_const(ind):
        return sum([isfloat(node) for node in ind])

    n_pop = len(pop)

    result_dict = dict()
    pop = [[node.name for node in ind]for ind in pop]
    
    for d in func_name:
        result_dict[d] = len([1 for ind in pop if d in ind])/n_pop
    for e, d in enumerate(x_name):
        result_dict[d] = len([1 for ind in pop if f'ARG{e}' in ind])/n_pop
    c_list = [count_const(ind) for ind in pop]    
    result_dict['const'] = sum([bool(i) for i in c_list])/n_pop
    result_dict['const/expr'] = sum(c_list)/n_pop
    return result_dict


