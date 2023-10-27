import decimal
import itertools
import operator
import os
import pickle
import random
import re
import time
import warnings
import copy
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# import pydotplus
import sympy
from deap import algorithms, base, creator, gp, tools
import deap
# from IPython.display import Image
from scipy import optimize
from scipy.optimize import basinhopping, least_squares, leastsq, minimize
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import cross_val_predict

from .deap_based_FGP_algorithms import FGP_NLS_algorithm
from .deap_based_func_node import Node_space, NumpyBasedFunction
# , add, cube, div, exp, ln, mul, sqrt, square, sub, protected_division, protected_sqrt, protected_ln, 
from .log_manager import table_log, txt_log
from .my_plot import line_plot, scatter_plot

warnings.simplefilter('ignore', RuntimeWarning)

CAN_USE_FUNC = ['add', 'sub', 'mul', 'div', 'sqrt', 'square', 'cube', 'ln', 'exp', 'protected_division', 'protected_ln', 'protected_sqrt']

def isfloat(_str):
    try:
        float(_str)
    except:
        return False
    return True


def get_depth_list(individual):
    depth = list()
    depth_pool = [0]
    for _node in individual:
        current_depth = depth_pool.pop()
        depth_pool.extend([current_depth+1]*_node.arity)
        depth.append(current_depth)
    return depth


def FV_filter(individual, function_group, function_filter, variable_filter):
    if function_filter or variable_filter:
        depth = get_depth_list(individual)
        func_pool, vals_pool = list(), list()
        com_flat = list(itertools.chain.from_iterable(function_group))
        old_d = -1
        for e, now_d in enumerate(depth):
            _node = individual[e]
            _name = individual[e].name
            _arity = individual[e].arity
            
            if old_d >= now_d:
                for _ in range(abs(old_d - now_d)+1):
                    func_pool.pop()
                    
            if function_filter:
                if _name in com_flat:
                    com_bool = [_name in c for c in function_group]
                    consider_pairs_n = [c for b, c in zip(com_bool, function_group) if b][0]
                    if sum([(n in func_pool) for n in consider_pairs_n])!=0:
                        return False, '=>>F-error'

            if variable_filter: 
                if _arity == 0: # check x or const node
                    if isinstance(_node.value, float): # True -> const node
                        pass
                    else:
                        if _name in vals_pool:
                            return False, '=>>V-error'
                        else:
                            vals_pool.append(_name)

            func_pool.append(_name)
            old_d = now_d

        state = ''
        if variable_filter:
            state += '=>>V-pass-' + '-'.join(vals_pool)
        else:
            state += '=>>V-none'
        if function_filter:
            state += '=>>F-pass'
        else:
            state += '=>>F-none'
        return True, state

    else:
        return True, '=>>F-none=>>V-none'


from imblearn.over_sampling import SMOTE, SMOTENC

def smote(X, n):
    """
    Sampling based on Smote
    """
    x  = X.values
    sm = SMOTE(random_state=0, k_neighbors=5)
    y  = np.ones(len(x))
    y2 = np.zeros(n+len(x))
    x2 = np.zeros((n+len(x),x.shape[1]))
    xcat = np.vstack([x, x2])
    ycat = np.concatenate([y,y2])
    xnew,ynew =  sm.fit_resample(xcat, ycat)
    xaug = xnew[ynew==1]
    return xaug

    

def D_filter(x_domain, y_domain, y_pred, equal, xydomain_filter):
    if xydomain_filter:
        if (x_domain is None or y_domain is None or y_pred is None):
            raise NameError(f'When xydomain_filter = True, x_domain needs a dataframe and y_domain needs a tuple of minimum and maximum values.\n x_domain = {type(x_domain)}, y_domain = {y_domain}')
        
        y_domain_min, y_domain_max = min(y_domain), max(y_domain)
        equal = ['=' if _e else '' for _e in equal]
        if eval(f'~np.all((y_domain_min <{equal[0]} y_pred)&(y_pred <{equal[0]} y_domain_max))'):
            return False, '=>>D-error'
        return True, '=>>D-pass'
    else:
        return True, '=>>D-none'

def C_only_filter(individual, constonly_filter):
    if constonly_filter:
        terminals=[isfloat(node.name) for node in individual if node.arity==0]
        if sum(terminals) == len(terminals):
            return False, '=>>C-error'
        else:
            return True, '=>>C-pass'
    else:
        return True, '=>>C-none'

def optimize_equations(ind, X, Xnum, compiler):
    # compiler est.toolbox_.compile
    func = compiler(ind)
    x0 = X.iloc[Xnum].values

    def _func_min(x0):
        return func(*x0)

    def _func_max(x0):
        return -func(*x0)

    _debug = True
    def _dump_vars(_bounds, _x0):
        print("============ DUMP ============")
        print(_bounds)
        print(_x0)

    #bounds = [(mn, mx) for mn, mx in zip(X.min(), X.max())]
    bounds = [(mn - (mx - mn) *0.1, mx + (mx - mn) *0.1) for mn, mx in zip(X.min(), X.max())]
    try:
        result_min = minimize(_func_min, x0, method="L-BFGS-B", bounds=bounds)
    except ValueError as e:
        print("ValueError (optimize_min)", e)
        result_min = type("Hoge", (object,), {"fun": -np.inf})
        if _debug:
            _dump_vars(bounds, x0)
            raise ValueError("optimize_min")
            

    try:
        result_max = minimize(_func_max, x0, method="L-BFGS-B", bounds=bounds)
    except ValueError as e:
        print("ValueError (optimize_max)", e)
        result_max = type("Hoge", (object,), {"fun": -np.inf})
        if _debug:
            _dump_vars(bounds, x0)
            raise ValueError("optimize_max")
            
        
    return result_min.fun, -result_max.fun

def D_filter2(
    individual, 
    x_domain, 
    y_domain, 
    equal, 
    domain_filter,
    compiler
    ):

    if domain_filter:
        y_domain_min, y_domain_max = min(y_domain), max(y_domain)
        minmax_pool = []
        n_samples = x_domain.shape[0]
        x0_index_list = random.sample(range(n_samples), k=5)
        for Xnum in x0_index_list:
            minmax_pool.extend(
                optimize_equations(
                    individual, 
                    x_domain, 
                    Xnum, 
                    compiler
                    ))
        y_min = min(minmax_pool)
        y_max = max(minmax_pool)
        # print(individual, f' {y_min} {y_max}')
        if ((y_domain_min <= y_min) & (y_max <= y_domain_max)):
            return True, f'=>>D2-pass{(y_min, y_max)}-n{len(set(minmax_pool))}'
        return False, f'=>>D2-error{(y_min, y_max)}-n{len(set(minmax_pool))}'

    else:
        return True, '=>>D2-none'
    

def FVD_filter(
    individual, 
    function_filter = True, 
    variable_filter = True, 
    xydomain_filter = True, 
    constonly_filter = True,
    function_group=[ 
        ['sqrt', 'protected_sqrt'], 
        ['square', 'cube'], 
        ['ln', 'exp', 'protected_ln']],
    x_domain=None, 
    y_domain=None, 
    y_pred=None, 
    equal=(True, True),
    # x_train=None,
    domain_filter=False, 
    compiler=None    
    ):

    state = ''

    _bool, _state = C_only_filter(individual, constonly_filter)
    state += _state
    if _bool == False:
        return _bool, state
    
    _bool, _state = FV_filter(individual, function_group, function_filter, variable_filter)
    state += _state
    if _bool == False:
        return _bool, state

    _bool, _state = D_filter(x_domain, y_domain, y_pred, equal, xydomain_filter)
    state += _state
    if _bool == False:
        return _bool, state

    _bool, _state = D_filter2(individual, x_domain, y_domain, equal, domain_filter, compiler)
    state += _state
    if _bool == False:
        return _bool, state

    return True, state


class surveyed_individuals():
    def __init__(self, x_df):
        self.sympy_space_ = Node_space(x_df, func_list=CAN_USE_FUNC)
        self._surveyed_ind_pool_ = set([])
    
    def is_unobserved(self, ind, add_pool=True):
        expr = str(ind)
        if expr not in self._surveyed_ind_pool_:
            if add_pool: self._surveyed_ind_pool_.add(expr)
            return True
        else:
            return False


def gnoise(x):
    _std = np.std(x, axis=0)
    _ret = pd.DataFrame(index=x.index, columns=x.columns, dtype="float64")

    for row in x.index:
        _gnoise = np.random.normal(loc=0, scale=1, size=_ret.shape[1])
        _ret.loc[row, :] = pd.Series(_std*_gnoise, index=_ret.columns).astype("float64")
    
    # print(_ret)
    return _ret


class Symbolic_Reg(BaseEstimator, RegressorMixin):
    def __init__(self, 
                 population_size = 1000, 
                 generations     = 200, 
                 tournament_size = 5,
                 num_elite_select = 1,
                 max_depth       = 4,
                 function_set    = ['add', 'sub', 'mul', 'div', 'ln', 'exp', 'sqrt', 'square', 'cube'], 
                 metric          = 'mae',
                 p_crossover     = 0.7,
                 p_mutation      = 0.2,
                 random_state    = 0,
                 const_range     = (-1, 1),
                 init_max_trial  = 50000,
                 init_unique     = True,
                 var_max_trial   = 20,
                 function_filter = True, 
                 variable_filter = True, 
                 xydomain_filter = True,
                 constonly_filter= True,
                 domain_filter   = False,
                 dfilter_aug     = False,
                 x_domain        = None,
                 y_domain        = None,
                 domain_equal    = (True, True),
                 results_dir     = './deap_based_FGP_results',
                 stabilize       = 0,
                 s_gnoise        = False,
                 s_lmd1          = 1.0,
                 s_lmd2          = 0.5,
                 s_clmd1         = 1.0,
                 s_clmd2         = 0.1,
                 ):
        """[summary]

        Args:
            population_size (int, optional): [description]. Defaults to 1000.
            generations (int, optional): [description]. Defaults to 200.
            tournament_size (int, optional): [description]. Defaults to 5.
            num_elite_select (int, optional): [description]. Defaults to 1.
            max_depth (int, optional): [description]. Defaults to 4.
            function_set (list, optional): [description]. Defaults to ['add', 'sub', 'mul', 'div', 'ln', 'log', 'sqrt'].
            metric (str, optional): [description]. Defaults to 'mae'.
            p_crossover (float, optional): [description]. Defaults to 0.7.
            p_mutation (float, optional): [description]. Defaults to 0.2.
            random_state (int, optional): [description]. Defaults to 0.
            const_range (tuple, optional): [description]. Defaults to (0, 1).
            x_domain ([type], optional): [description]. Defaults to None.
            y_domain ([type], optional): [description]. Defaults to None.
            results_dir (str, optional): [description]. Defaults to './results'.
        """
        
        self.population_size = population_size
        self.generations = generations
        self.tournament_size = tournament_size
        self.num_elite_select = num_elite_select
        self.max_depth = max_depth
        self.function_set = function_set
        self.metric = metric
        self.p_crossover = p_crossover
        self.p_mutation = p_mutation
        self.random_state = random_state
        self.init_max_trial = init_max_trial
        self.init_unique = init_unique
        self.var_max_trial = var_max_trial
        self.const_range = const_range
        self.function_filter = function_filter
        self.variable_filter = variable_filter
        self.xydomain_filter = xydomain_filter
        self.constonly_filter= constonly_filter
        self.domain_filter   = domain_filter
        self.dfilter_aug     = dfilter_aug
        self.x_domain = x_domain
        self.y_domain = y_domain
        self.domain_equal = domain_equal
        
        self.results_dir = results_dir

        self._n_time = 1


        random.seed(self.random_state)
        self.fit_x_ = None
        self.fit_y_ = None
        os.makedirs(self.results_dir, exist_ok=True)
        self.text_log_ = txt_log(file_name='000_GP_log_txt', save_path=self.results_dir)
        
        self._n_ind_generations = 0
        self._n_ind_gen_successes = 0
        self._warnings = 0
        
        self._can_use_func = ['add', 'sub', 'mul', 'div', 'sqrt', 'square', 'cube', 'ln', 'exp', 'protected_division', 'protected_ln', 'protected_sqrt']

        """stabilize"""
        self.stabilize = stabilize
        self.s_gnoise = s_gnoise
        self.s_lmd1 = s_lmd1
        self.s_lmd2 = s_lmd2
        self.s_clmd1 = s_clmd1
        self.s_clmd2 = s_clmd2


    def fit(self, x, y):
        _start_gen = time.time()
        self.fit_x_ = x
        self.fit_y_ = y
        
        self._surveyed_individuals_ = surveyed_individuals(self.fit_x_)
        
        self.pset = gp.PrimitiveSet("MAIN", x.shape[1])
        
        for i, x_name in enumerate(x.columns):
            p = {'ARG{}'.format(i):f'{x_name}'}
            self.pset.renameArguments(**p)

        if 'add' in self.function_set: self.pset.addPrimitive(NumpyBasedFunction.add, 2)
        if 'sub' in self.function_set: self.pset.addPrimitive(NumpyBasedFunction.sub, 2)
        if 'mul' in self.function_set: self.pset.addPrimitive(NumpyBasedFunction.mul, 2)
        if 'div' in self.function_set: self.pset.addPrimitive(NumpyBasedFunction.div, 2)
        if 'ln'  in self.function_set: self.pset.addPrimitive(NumpyBasedFunction.ln, 1)
        if 'sqrt' in self.function_set: self.pset.addPrimitive(NumpyBasedFunction.sqrt, 1)
        if 'square' in self.function_set: self.pset.addPrimitive(NumpyBasedFunction.square, 1)
        if 'cube' in self.function_set: self.pset.addPrimitive(NumpyBasedFunction.cube, 1)
        if 'exp' in self.function_set: self.pset.addPrimitive(NumpyBasedFunction.exp, 1)
        if 'protected_division' in self.function_set: self.pset.addPrimitive(NumpyBasedFunction.protected_division, 2)
        if 'protected_sqrt' in self.function_set: self.pset.addPrimitive(NumpyBasedFunction.protected_sqrt, 1)
        if 'protected_ln' in self.function_set: self.pset.addPrimitive(NumpyBasedFunction.protected_ln, 1)

    
        # add initial constant to be optimized
        n_c_node = 1
        add_n_c_node = 0
        _count = 0
        while add_n_c_node < n_c_node:
            try:
                self.pset.addEphemeralConstant(f'c_node_{_count}', lambda: random.uniform(self.const_range[0],self.const_range[1]))
                add_n_c_node += 1
                _count += 1
            except:
                _count += 1
                pass

        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
        creator.create("Individual", gp.PrimitiveTree, fitness=creator.FitnessMin)

        def _progress_bar():
            if int(self.population_size*0.01)*self._n_time == self._n_ind_gen_successes:
                if self._n_time%10 == 0:
                    print(f'{int(self._n_time)}%', end='')
                else:
                    print('|', end='')
                self._n_time += 1

        def filter_initIterate(container, generator, max_trial=self.init_max_trial, unique=self.init_unique, text_log=self.text_log_):
            for i in range(max_trial):
                ind = creator.Individual(generator())
                score = self._evalSymbReg(ind, self.fit_x_, self.fit_y_)
                self._n_ind_generations += 1
                
                if score[0] != np.inf:
                    if unique:
                        if self._surveyed_individuals_.is_unobserved(creator.Individual(ind)):
                            self._n_ind_gen_successes += 1
                            _progress_bar()
                            return container(ind)
                    else:
                        self._n_ind_gen_successes += 1
                        _progress_bar()
                        return container(ind)
                
            raise NameError(f'The maximum number of trials has been reached. \nNumber already generated : {self._n_ind_gen_successes}\nNumber of challenges : {self._n_ind_generations}')

        self.toolbox_ = base.Toolbox()
        self.toolbox_.register("expr", gp.genHalfAndHalf, pset=self.pset, min_=1, max_=2) # gp.genHalfAndHalf https://deap.readthedocs.io/en/master/api/tools.html#deap.gp.genHalfAndHalf
        
        self.toolbox_.register("individual", filter_initIterate, creator.Individual, self.toolbox_.expr)
        # if (self.function_filter | self.variable_filter | self.xydomain_filter):
        #     self.toolbox_.register("individual", filter_initIterate, creator.Individual, self.toolbox_.expr)
        # else:
        #     # Normal generation without filter
        #     self.toolbox_.register("individual", tools.initIterate, creator.Individual, self.toolbox_.expr)
            
        self.toolbox_.register("population", tools.initRepeat, list, self.toolbox_.individual)
        self.toolbox_.register("compile", gp.compile, pset=self.pset)
        self.toolbox_.register("evaluate", self._evalSymbReg, x=x, y_true=y)
        self.toolbox_.register("select", tools.selTournament, tournsize=self.tournament_size)

        # gp.cxOnePoint : 1 point crossover
        # https://deap.readthedocs.io/en/master/api/tools.html?highlight=bloat#deap.gp.cxOnePoint
        self.toolbox_.register("mate", gp.cxOnePoint) 
        self.toolbox_.decorate("mate", gp.staticLimit(key=operator.attrgetter("height"), max_value=self.max_depth))

        self.toolbox_.register("expr_mut", gp.genFull, min_=0, max_=2)

        # gp.mutUniform 
        # https://deap.readthedocs.io/en/master/api/tools.html?highlight=bloat#deap.gp.mutUniform
        self.toolbox_.register("mutate", gp.mutUniform, expr=self.toolbox_.expr_mut, pset=self.pset) 
        self.toolbox_.decorate("mutate", gp.staticLimit(key=operator.attrgetter("height"), max_value=self.max_depth))
        
        self._time0 = time.time()

        print('Generation of initial generation')
        pop = self.toolbox_.population(n=self.population_size)
        print('\nGeneration complete')
        
        self.text_log_.print([f' {self._n_ind_gen_successes} [ ind ] / {self._n_ind_generations} [ trials ] (time : {(time.time() - self._time0)/60:.3f} min)\n'])
        
        self.hof = tools.HallOfFame(self.num_elite_select)
        
        pop, log = FGP_NLS_algorithm(population       = pop,
                            toolbox          = self.toolbox_, 
                            cxpb             = self.p_crossover,
                            mutpb            = self.p_mutation,
                            ngen             = self.generations,
                            halloffame       = self.hof,
                            num_elite_select = self.num_elite_select,
                            var_max_trial    = self.var_max_trial, 
                            check_func       = self._surveyed_individuals_,
                            text_log         = self.text_log_,
                            save_dir         = self.results_dir,
                            func_name        = self.function_set)
        
        self.expr = tools.selBest(pop, 1)[0]
        self.tree = gp.PrimitiveTree(self.expr)
        self.nodes, self.edges, self.labels = gp.graph(self.expr)
        self.log = log
        
        self.text_log_.print([f'FGP-NLS All Execution Time : {time.time() - _start_gen:.3f} s', 
                              f'FGP-NLS All Execution Time : {(time.time() - _start_gen)/60:.1f} min',
                              f'FGP-NLS All Execution Time : {(time.time() - _start_gen)/60/60:.1f} h',
                              f'Number of constant optimization warnings : {self._warnings}'
                              ])
        return self

    def predict(self, x):
        y_pred = self._pred(x, self.expr)
        return y_pred
    
    def _pred(self, x, expr):
        func = self.toolbox_.compile(expr=expr)
        x_data = (x['{}'.format(i)] for i in list(x.columns))
        # y_pred = func(*x_data)
        try:
            y_pred = func(*x_data)
        except:
            self._warnings += 1
            self.text_log_.print(['ERROR!! nan is included.', f'{str(expr)}', self.root])
            self.text_log_.print(self.temporary)
            self.text_log_.print(['ERROR!! _results'])
            self.text_log_.print(self.temporary2)
            
            y_pred = np.inf 
            
        if np.isscalar(y_pred):      # Avoid scalar errors.
            y_pred = pd.Series(np.full(x.shape[0], float(y_pred)))
        elif len(y_pred.shape) == 0: # Avoid errors due to singleton arrays.
            y_pred = y_pred.item()
            y_pred = pd.Series(np.full(x.shape[0], float(y_pred)))
        return y_pred
        
    def _evalSymbReg(self, individual, x, y_true):
        individual.state = ''

        # >>>>> func of const opt
        def _func(constants, x, y, individual, compiler, constant_nodes):
            _idx = 0
            for i in constant_nodes:
                if ~np.isnan(constants[_idx]):
                    cnode = copy.deepcopy(individual[i])
                    cnode.value = constants[_idx]
                    cnode.name  = str(constants[_idx])
                    
                    individual[i] = cnode
                _idx += 1

            _f = compiler(expr=individual)
            _x_data = (x['{}'.format(i)] for i in list(x.columns))
            _y = _f(*_x_data)
            if np.isscalar(_y):         # Avoid scalar errors.
                _y = np.full(x.shape[0], float(_y))
            elif len(_y.shape) == 0:    # Avoid errors due to singleton arrays.
                _y = _y.item()
                _y = np.full(x.shape[0], float(_y))
            elif type(_y) is pd.Series: # Set it to np.ndarray
                _y = _y.values
            y = y.values.reshape(len(_y))
            residual = y - _y
            
            return residual
        # <<<<< func of const opt
        
        filter_results1 = FVD_filter(individual, 
                                    function_filter = self.function_filter, 
                                    variable_filter = self.variable_filter, 
                                    xydomain_filter = False, 
                                    constonly_filter = self.constonly_filter,
                                    function_group=[ 
                                        ['sqrt'], 
                                        ['square', 'cube'], 
                                        ['ln', 'log', 'exp']],
                                    x_domain=None, 
                                    y_domain=None, 
                                    y_pred=None, 
                                    equal=self.domain_equal,
                                    # x_train=None,
                                    domain_filter=False, 
                                    compiler=None                                    
                                    )
        if filter_results1[0]:
            pass
        else:
            individual.state = filter_results1[1]
            return np.inf,

        _is_const = [isfloat(n.name) for n in individual]

        if sum(_is_const):
            constant_nodes = [e for e, i in enumerate(_is_const) if i]
            constants0 = [individual[idx].value for idx in constant_nodes]
            self.temporary = constants0
            self.temporary2 = ['']
            if 0 < len(constants0):
                try:
                    _result=least_squares(_func, x0 = constants0, args=(x, y_true, individual, self.toolbox_.compile, constant_nodes), method='lm')
                    _idx = 0
                    if _result.status >= 2:
                        for i in constant_nodes:
                            cnode = copy.deepcopy(individual[i])
                            cnode.value = _result.x[_idx]
                            cnode.name  = str(_result.x[_idx])
                            individual[i] = cnode
                            _idx += 1
                        self.root = 'A'
                        self.temporary2 = [_result.x, _result.success, 'status', _result.status, _result.message]
                        
                    else:
                        for i in constant_nodes:
                            cnode = copy.deepcopy(individual[i])
                            if ~np.isnan(constants0[_idx]):
                                cnode.value = constants0[_idx]
                                cnode.name  = str(constants0[_idx])
                            
                            individual[i] = cnode
                            _idx += 1
                        self.root = 'B'
                        
                except:
                    _idx = 0
                    for i in constant_nodes:
                        cnode = copy.deepcopy(individual[i])
                        if ~np.isnan(constants0[_idx]):
                            cnode.value = constants0[_idx]
                            cnode.name  = str(constants0[_idx])
                            
                        individual[i] = cnode
                        _idx += 1
                    self.root = 'C'

        if self.dfilter_aug:
            _x_domain = pd.DataFrame(smote(self.x_domain, 10000), columns=self.x_domain.columns)
            _dfilter_aug_log = '=>>X-Aug'
        else:
            _x_domain = self.x_domain
            _dfilter_aug_log = ''
            
        #print("smote")
        filter_results2 = FVD_filter(individual, 
                                    function_filter = False, 
                                    variable_filter = False, 
                                    xydomain_filter = self.xydomain_filter,
                                    constonly_filter = False,
                                    # x_domain=self.x_domain, 
                                    x_domain=_x_domain,
                                    y_domain=self.y_domain, 
                                    # y_pred=self._pred(self.x_domain, individual), 
                                    y_pred=self._pred(_x_domain, individual), 
                                    equal=self.domain_equal,
                                    # x_train=x,
                                    domain_filter= self.domain_filter, 
                                    compiler=self.toolbox_.compile
                                    )
        if filter_results2[0]:
            pass
        else:
            individual.state = filter_results1[1] + _dfilter_aug_log + filter_results2[1]
            return np.inf,
        
        y_pred = self._pred(self.fit_x_, individual)
        individual.state = filter_results1[1] + _dfilter_aug_log + filter_results2[1]

        try:
            if self.metric == 'mae':
                error = mean_absolute_error(y_true, y_pred)
            elif self.metric == 'rmse':
                error = mean_squared_error(y_true, y_pred, squared=False)
            elif self.metric == 'mse':
                error = mean_squared_error(y_true, y_pred, squared=True)
            else:
                error = mean_absolute_error(y_true, y_pred)
        except:
            individual.state += '=opt-error'
            error = np.inf
        
        """CHECK STABILITY
        """

        """stability for variable noise"""
        if 1 == self.stabilize or 3 == self.stabilize:
            error_noise = None
            
            if self.s_gnoise:
                # _x_gnoise = self.fit_x_
                _x_gnoise = self.fit_x_ + self.s_lmd2*gnoise(self.fit_x_)
                # print(self.fit_x_)
                # print(_x_gnoise)
                y_pred_noise = self._pred(_x_gnoise, individual)
                try:
                    if self.metric == 'mae':
                        error_noise = mean_absolute_error(y_pred, y_pred_noise)
                    elif self.metric == 'rmse':
                        error_noise = mean_squared_error(y_pred, y_pred_noise, squared=False)
                    elif self.metric == 'mse':
                        error_noise = mean_squared_error(y_pred, y_pred_noise, squared=True)
                    else:
                        error_noise = mean_absolute_error(y_pred, y_pred_noise)
                except:
                    error_noise = np.inf
            else:
                # _x_gnoise = self.fit_x_
                _std = np.std(self.fit_x_, axis=0)
                _x_noise_p = self.fit_x_ + self.s_lmd2*_std
                _x_noise_m = self.fit_x_ - self.s_lmd2*_std
                # print("org")
                # print(self.fit_x_.iloc[0,:])
                # print("noise p")
                # print(_x_noise_p.iloc[0,:])
                # print("noise p")
                # print(_x_noise_m.iloc[0,:])
                y_pred_noise = self._pred(pd.concat([_x_noise_p,_x_noise_m]), individual)
                # print(pd.concat([_x_noise_p,_x_noise_m]).shape, _x_noise_p.shape)

                try:
                    if self.metric == 'mae':
                        error_noise = mean_absolute_error(pd.concat([y_pred,y_pred]), y_pred_noise)
                    elif self.metric == 'rmse':
                        error_noise = mean_squared_error(pd.concat([y_pred,y_pred]), y_pred_noise, squared=False)
                    elif self.metric == 'mse':
                        error_noise = mean_squared_error(pd.concat([y_pred,y_pred]), y_pred_noise, squared=True)
                    else:
                        error_noise = mean_absolute_error(pd.concat([y_pred,y_pred]), y_pred_noise)
                except:
                    error_noise = np.inf

            if error_noise is not None and self.s_lmd1 != 0:
                error += self.s_lmd1*error_noise
                # return error + self.stb_lambda*error_noise, error, error_noise


        """stability for regression coefficient"""
        if 2 == self.stabilize or 3 == self.stabilize:
            error_noise = None

            #print(individual)
            _org_coeffs = {}
            for idx,node in enumerate(individual):
                if node.arity == 0 and isfloat(node.name):
                    #print(individual, individual[idx].value)
                    _org_coeffs[idx] = node.value
                    #print(type(individual[idx].value))
                    individual[idx].value = node.value*(1.0 + self.s_clmd2)
                    #print(individual, individual[idx].value)
                    #print(_org, individual[idx].value, individual[idx].value-_org)

            #print(individual)
            y_pred_noise_p = self._pred(self.fit_x_, individual)
            """restore org coeffs"""
            for idx,value in _org_coeffs.items():
                individual[idx].value = value

            _org_coeffs = {}
            for idx,node in enumerate(individual):
                if node.arity == 0 and isfloat(node.name):
                    _org = node.value
                    _org_coeffs[idx] = node.value
                    individual[idx].value = node.value*(1.0 - self.s_clmd2)
                    #print(_org, individual[idx].value, individual[idx].value-float(_org))
                    
            #print(individual)
            y_pred_noise_m = self._pred(self.fit_x_, individual)
            """restore org coeffs"""
            for idx,value in _org_coeffs.items():
                individual[idx].value = value
                
            #print(y_pred.iloc[0], y_pred_noise_p.iloc[0], y_pred_noise_m.iloc[0])
                

            try:
                if self.metric == 'mae':
                    error_noise = mean_absolute_error(pd.concat([y_pred,y_pred]), pd.concat([y_pred_noise_p, y_pred_noise_m]))
                elif self.metric == 'rmse':
                    error_noise = mean_squared_error(pd.concat([y_pred,y_pred]), pd.concat([y_pred_noise_p, y_pred_noise_m]), squared=False)
                elif self.metric == 'mse':
                    error_noise = mean_squared_error(pd.concat([y_pred,y_pred]), pd.concat([y_pred_noise_p, y_pred_noise_m]), squared=True)
                else:
                    error_noise = mean_absolute_error(pd.concat([y_pred,y_pred]), pd.concat([y_pred_noise_p, y_pred_noise_m]))
            except:
                error_noise = np.inf

            if error_noise is not None and self.s_clmd1 != 0:
                error += self.s_clmd1*error_noise
                # return error + self.stb_lambda*error_noise, error, error_noise

            
        return error,
    
    def save_gen_metric_plot(self):
        log_file = pd.DataFrame(self.log)
        log_file.to_csv(f'{self.results_dir}/001_GP_log.tsv', sep='\t')
        
        line_plot(x_data=list(log_file.index), y_data=log_file.loc[:, ['score-min']], 
                  c_data=['b', 'r', 'g'], label_data=[f'{self.metric}-min'], 
                  xy_labels=['generation', f'{self.metric}'], figsize=(16,8), 
                  cmap='jet', color_bar=True, save_name=f'{self.results_dir}/001_GP_log_min', data_direction='column',
                  invert_xaxis=False, linewidth = 2.0, font_size=30, facecolor=None, vspan_data=None, line_style=['-', '--', ':'], )
        
        line_plot(x_data=list(log_file.index), y_data=log_file.loc[:, [f'score-min', f'score-med']], c_data=['b', 'r', 'g'], label_data=[f'{self.metric}-min', f'{self.metric}-med'], xy_labels=['generation', f'{self.metric}'], figsize=(16,8), 
                    cmap='jet', color_bar=True, save_name=f'{self.results_dir}/001_GP_log', data_direction='column', 
                    invert_xaxis=False, linewidth = 2.0, font_size=30, facecolor=None, vspan_data=None, line_style=['-', '--', ':'])
        line_plot(x_data=list(log_file.index), y_data=log_file.loc[:, [f'score-std']], c_data=['g'], label_data=[f'{self.metric}-std'], xy_labels=['generation', f'{self.metric}'], figsize=(16,8), 
                    cmap='jet', color_bar=True, save_name=f'{self.results_dir}/001_GP_log_std', data_direction='column', 
                    invert_xaxis=False, linewidth = 2.0, font_size=30, facecolor=None, vspan_data=None, line_style=['-', '--', ':'])
        line_plot(x_data=list(log_file.index), y_data=log_file.loc[:, ['unique_rate']], c_data=['m'], label_data=['unique_rate'], xy_labels=['generation', 'unique_rate'], figsize=(16,8), 
                    cmap='jet', color_bar=True, save_name=f'{self.results_dir}/001_GP_log_unique_rate', data_direction='column', 
                    invert_xaxis=False, linewidth = 2.0, font_size=30, facecolor=None, vspan_data=None, line_style=['-', '--', ':'])
        return log_file

    def save_expr4word(self):
        sympy_space = Node_space(self.fit_x_, func_list=self._can_use_func)   
        expr = eval(str(self.expr), sympy_space.Nspace())
        expr4word = sympy.printing.octave.octave_code(expr).replace('./', '/').replace('.^', '^').replace('.*', '\\cdot ').replace('*', ' \\cdot ').replace('sqrt', '\\sqrt')
        f = open(f'{self.results_dir}/002_best_expr4word.txt', 'w', encoding='cp932')
        f.write(expr4word)
        f.close()
        
        non_e_expr4word = expr4word.replace('e+', '\\cdot 10^').replace('e-', '\\cdot 10^-')
        f = open(f'{self.results_dir}/002_best_expr4word_non_e.txt', 'w', encoding='cp932')
        f.write(non_e_expr4word)
        f.close()
        
        num_only = re.findall(r"\d+\.\d+", expr4word)
        for num in num_only:
            expr4word = expr4word.replace(num, Significant_figures(num, digit=3))
        
        f = open(f'{self.results_dir}/002_best_expr4word_3digits.txt', 'w', encoding='cp932')
        f.write(expr4word)
        f.close()
        
        non_e_expr4word = expr4word.replace('e+', '\\cdot 10^').replace('e-', '\\cdot 10^-')
        f = open(f'{self.results_dir}/002_best_expr4word_3digits_non_e.txt', 'w', encoding='cp932')
        f.write(non_e_expr4word)
        f.close()
        return expr4word
    
    def save_expr4tex(self, y_name=''):
        sympy_space = Node_space(self.fit_x_, func_list=self._can_use_func)
        # sympy.var(self.node_space_.symbol)
        expr = eval(str(self.expr), sympy_space.Nspace())
        expr_ = sympy.latex(expr)
        expr_la = expr_.replace('\\\\', '\\')
        f = open(f'{self.results_dir}/002_best_model.tex', 'w', encoding='cp932')
        f.write('\documentclass[a4paper, 12pt, fleqn]{tarticle}\n')
        f.write('\\begin{document}\n')
        f.write('\\begin{equation}\n')
        f.write(y_name + expr_la + '\n')
        f.write('\\end{equation}\n')
        f.write('\\end{document}\n')
        f.close()

        expr4png = sympy.latex(expr, mul_symbol='dot')#, min=-4, max=4)
        # expr4png = sympy.printing.latex(expr, mul_symbol='dot', min=-4, max=4)
        decimal_list = re.findall(r"\d\.\d+", expr4png)
        for string in decimal_list:
            if -1 < float(string) < 1:
                expr4png = expr4png.replace(string, f'{float(string):.4f}', 1)
            elif (-10 < float(string) <= -1 or 1 <= float(string) < 10):
                expr4png = expr4png.replace(string, f'{float(string):.3f}', 1)
            elif (-100 < float(string) <= -10 or 10 <= float(string) < 100):
                expr4png = expr4png.replace(string, f'{float(string):.2f}', 1)
            elif (-1000 < float(string) <= -100 or 100 <= float(string) < 1000):
                expr4png = expr4png.replace(string, f'{float(string):.1f}', 1)
            else:
                expr4png = expr4png.replace(string, f'{float(string):.0f}', 1)

        expr4png = y_name + expr4png
        tex = f'${expr4png}$'

        fig, ax = plt.subplots(figsize=(1, 1))
        ax.text(0.4, 0.4, rf'{tex}', size=30)
        plt.axis("off")
        plt.savefig(f'{self.results_dir}/002_best_model_expr.png', dpi=400, transparent=True, bbox_inches="tight", pad_inches=0.05)
        plt.close()
        return expr4png
    
    def save_node_analysis(self):
        log_file = pd.read_csv(f'{self.results_dir}/001_GP_node_analysis.tsv', sep='\t', index_col=0)
        _y_data = log_file.loc[:, list(self.fit_x_.columns)]
        line_plot(x_data=list(log_file.index), y_data=_y_data, c_data=list(range(len(_y_data.columns))), label_data=list(_y_data.columns), xy_labels=['generation', 'Content rate'], figsize=(16,8), 
                    color_bar=False, save_name=f'{self.results_dir}/001_GP_node_analysis_X', data_direction='column', 
                    invert_xaxis=False, linewidth = 2.0, font_size=30, cmap='my_cmap', line_style=['-', '--', ':'])
        log_file = log_file.drop(list(self.fit_x_.columns), axis=1)
        line_plot(x_data=list(log_file.index), y_data=log_file, c_data=list(range(len(log_file.columns))), label_data=list(log_file.columns), xy_labels=['generation', 'Content rate'], figsize=(16,8), 
                    color_bar=False, save_name=f'{self.results_dir}/001_GP_node_analysis_func_all', data_direction='column', 
                    invert_xaxis=False, linewidth = 2.0, font_size=30, cmap='my_cmap', line_style=['-', '--', ':'])
        log_file = log_file.drop('const/expr', axis=1)
        line_plot(x_data=list(log_file.index), y_data=log_file, c_data=list(range(len(log_file.columns))), label_data=list(log_file.columns), xy_labels=['generation', 'Content rate'], figsize=(16,8), 
                    color_bar=False, save_name=f'{self.results_dir}/001_GP_node_analysis_func', data_direction='column', 
                    invert_xaxis=False, linewidth = 2.0, font_size=30, cmap='my_cmap', line_style=['-', '--', ':'])
        
        log_file = pd.read_csv(f'{self.results_dir}/001_GP_node_analysis_select.tsv', sep='\t', index_col=0)
        _y_data = log_file.loc[:, list(self.fit_x_.columns)]
        line_plot(x_data=list(log_file.index), y_data=_y_data, c_data=list(range(len(_y_data.columns))), label_data=list(_y_data.columns), xy_labels=['generation', 'Content rate'], figsize=(16,8), 
                    color_bar=False, save_name=f'{self.results_dir}/001_GP_node_analysis_X_select', data_direction='column', 
                    invert_xaxis=False, linewidth = 2.0, font_size=30, cmap='my_cmap', line_style=['-', '--', ':'])
        log_file = log_file.drop(list(self.fit_x_.columns), axis=1)
        line_plot(x_data=list(log_file.index), y_data=log_file, c_data=list(range(len(log_file.columns))), label_data=list(log_file.columns), xy_labels=['generation', 'Content rate'], figsize=(16,8), 
                    color_bar=False, save_name=f'{self.results_dir}/001_GP_node_analysis_func_all_select', data_direction='column', 
                    invert_xaxis=False, linewidth = 2.0, font_size=30, cmap='my_cmap', line_style=['-', '--', ':'])
        log_file = log_file.drop('const/expr', axis=1)
        line_plot(x_data=list(log_file.index), y_data=log_file, c_data=list(range(len(log_file.columns))), label_data=list(log_file.columns), xy_labels=['generation', 'Content rate'], figsize=(16,8), 
                    color_bar=False, save_name=f'{self.results_dir}/001_GP_node_analysis_func_select', data_direction='column', 
                    invert_xaxis=False, linewidth = 2.0, font_size=30, cmap='my_cmap', line_style=['-', '--', ':'])
    
    def save_gp_params(self):
        p_df = pd.Series(self.get_params())
        p_df.to_csv(f'{self.results_dir}/001_FilterGPSR_params.tsv', '\t')

    def save_expr(self):
        save_dict = dict(expr=str(self.expr), use_X=[node.value for node in self.expr if type(node)==deap.gp.Terminal])
        with open(f'{self.results_dir}/000_best_expr.json', 'w') as fp:
            json.dump(save_dict, fp)
    
    def save_all(self, y_name=''):
        self.save_gp_params()
        self.save_gen_metric_plot()
        self.save_expr4tex(y_name=y_name)
        self.save_expr4word()
        self.save_expr()
        self.save_node_analysis()
        # self.save_tree_pic()
        

def Significant_figures(num, digit=3):
    decimal.getcontext().prec = digit
    a = decimal.Decimal(num)
    b = decimal.Decimal(1)
    return str(a/b)

def output_score(   y_true_list,
                    y_pred_list,
                    data_name_list = ['train', 'test'],
                    save_name      = './'
                    ):
    metric = dict( R2 = r2_score, MAE = mean_absolute_error, RMSE = mean_squared_error)
    colname = []
    score = pd.DataFrame()
    for idx, _d in enumerate(y_true_list):
        data_name = data_name_list[idx]
        for key in metric:
            d_key = '{}_{}'.format(key, data_name)
            if key == 'RMSE':
                try:
                    score.at[0,d_key] = metric[key](y_true_list[idx], y_pred_list[idx], squared = False)
                except:
                    score.at[0,d_key] = np.inf
            elif key == 'R2':
                try:
                    score.at[0,d_key] = metric[key](y_true_list[idx], y_pred_list[idx])
                except:
                    score.at[0,d_key] = -np.inf
            else:
                try:
                    score.at[0,d_key] = metric[key](y_true_list[idx], y_pred_list[idx])
                except:
                    score.at[0,d_key] = np.inf
    score.to_csv('{}score.tsv'.format(save_name), sep='\t')
    


class ExprNameSpace():
    def __init__(self, X, func=['add', 'sub', 'mul', 'div', 'ln', 'sqrt', 'square', 'cube', 'exp']):
        self.X = X
        self.func = func
        self.nspace = {}
        self._make()
    
    def _make(self):
        self.nspace.update({_f:eval(f'NumpyBasedFunction.{_f}')for _f in self.func})
        self.nspace.update({xname:np.array(self.X[xname]).reshape(-1, 1) for xname in self.X.columns})


class LoadExpr():
    def __init__(self, expr):
        self.expr = expr

    def predict(self, X):
        ns = ExprNameSpace(X)
        return eval(str(self.expr), ns.nspace)

