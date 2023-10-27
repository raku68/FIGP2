import copy
import itertools
import random
import time
from collections import deque

from joblib import Parallel, delayed

import numpy as np
import pandas as pd
import sympy
from scipy.optimize import basinhopping, least_squares, leastsq, minimize
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


class tree_graph(list):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.optimized = False
        self._const_list, self._depth, self.score = None, None, None
        self._expr = ''
    
    def __setitem__(self, key, value):
        super().__setitem__(key, value)
        self.optimized = False
        self._const_list, self._depth, self.score = None, None, None
        self._expr = ''
        return self
    
    def __deepcopy__(self, memo):
        new = self.__class__(self)
        new.__dict__.update(copy.deepcopy(self.__dict__, memo))
        return new
    
    def append(self, object):
        super().append(object)
        self.optimized = False
        self._const_list, self._depth, self.score = None, None, None
        self._expr = ''
    
    def extend(self, iterable):
        super().extend(iterable)
        self.optimized = False
        self._const_list, self._depth, self.score = None, None, None
        self._expr = ''

    def calc_depth(self):
        if self._depth is None:
            _const_list, _depth = list(), list()
            depth_pool = [0]
            n_const = 0
            for e, _t in enumerate(self):
                current_depth = depth_pool.pop()
                depth_pool.extend([current_depth+1]*_t['arity'])
                _depth.append(current_depth)
                if _t['value'] is not None:
                    val = copy.deepcopy(_t['value'])
                    _const_list.append(val)
                    self[e] = {'name':f'const{n_const}', 'arity':0, 'value':val}
                    n_const += 1
            self._const_list = _const_list
            self._depth = _depth
        return self
    
    def update_value(self, const_list):
        old_depth = self._depth
        queue = deque(const_list)
        n_const = 0
        if len(self) > 0:
            for e, n in enumerate(self):
                if n['value'] is not None:
                    self[e] = {'name':f'const{n_const}', 'arity':0, 'value':queue.popleft()}
                    n_const += 1
        self.optimized = True
        self._const_list = const_list
        self._depth = old_depth
        return self
    
    def get_subtree(self, begin):
        end = begin + 1
        total = self[begin]['arity']
        while total > 0:
            total += self[end]['arity'] - 1
            end += 1
        return slice(begin, end)
    
    def to_str(self):
        if self._depth is None:
            self.calc_depth()
        self._expr = ''
        last = len(self)
        step = 0
        old_dept = -1
        count = 0
        for node, dept in zip(self, self._depth):#node, dept
            step += 1
            if old_dept == -1:
                self._expr += ''
            elif old_dept < dept:
                self._expr += '('
                count += 1
            elif old_dept == dept:
                self._expr += ','
            elif old_dept > dept:
                self._expr += ')' * abs(old_dept - dept) + ','
                count -= abs(old_dept - dept)
            self._expr += str(node['name'])
            old_dept = dept
            if step == last:
                self._expr += ')' * count
        return self._expr
    
    def pred(self, node_space):
        return eval(self.expr, node_space)

    def evaluation(self, node_space, y_true, X=None, metric='MAE', save=True, forced=False):
        if self.score is not None:
            if forced:
                pass
            else:
                return self.score
        y_true = np.array(y_true).reshape(-1)
        y_pred = self.pred(node_space.update(X=X, const_list=self._const_list))
        if isinstance(y_pred, np.float64):
            y_pred = np.full(y_true.shape, y_pred)
            
        if metric == 'MAE':
            try:
                score = mean_absolute_error(y_true, y_pred)
            except:
                score = np.inf

        elif metric == 'RMSE':
            try:
                score = mean_squared_error(y_true, y_pred, squared=False)
            except:
                score = np.inf

        elif metric == 'MSE':
            try:
                score = mean_squared_error(y_true, y_pred, squared=True)
            except:
                score = np.inf

        elif metric == 'R2':
            try:
                score = r2_score(y_true, y_pred)
            except:
                score = np.NINF
        if save:
            self.score = score
        return score

    def save_node_spase(self, node_space):
        self.node_spase_ = copy.deepcopy(node_space)

    def predict(self, X):
        return self.pred(self.node_spase_.update(X=X, const_list=self._const_list, use_numpy=True))


    @property
    def expr(self):
        if len(self._expr) == 0:
            self.to_str()
        return self._expr

    @property
    def const_list(self):
        return self._const_list
    



class make_tree():
    def __init__(self, node_space):
        self.node_space = node_space
    
    def pop_generate(   self, pop_size, depth, 
                        full=1, 
                        grow=1,
                        c_range=(-1,1), 
                        c_rate=None, 
                        c_opt = None,
                        y_true = None,
                        grow_50_50=False,
                        linear_scaling=False,
                        function_filter=False,
                        variable_filter=False, 
                        xydomain_filter=False,
                        fcom = [['sqrt'], ['Square', 'Cube'], ['log'], ['exp', 'Pow', 'Ipow']],
                        x_domain=None, 
                        y_domain=None, 
                        equal=(True, True),
                        txlog=None,
                        unique_only=False):

        base_time = time.time()
        if c_opt is not None:
            if y_true is None:
                raise Exception('Y_true is required for constant node optimization.')

        times = 2000

        txt_list = list()
        node_space = self.node_space
        domain_node_space_ = copy.deepcopy(self.node_space)

        if c_rate is None:
            c_rate = 1/(len(node_space.Xinfo) + 1)
        else:
            if isinstance(c_rate, float):
                pass
            else:
                raise Exception(f'c_rate must be None or float. {c_rate} was given.')

        n_full = int(full/(full+grow) * pop_size)
        n_grow = pop_size - n_full
        L1 = [i['name'] for i in node_space.Finfo]
        L2 = [i['name'] for i in node_space.Xinfo]
        txt_list.extend(['---- pop_generate start ----', '-- Parameters ', f'Use function_filter : {L1}', f'Number of functions : {len(node_space.Finfo)}', f'Use descriptor : {L2}', f'Number of descriptors : {len(node_space.Xinfo)}', f'grow_50_50 : {grow_50_50}', f'depth : {depth}'])
        txt_list.extend([f'c_range : {c_range}', f'c_rate : {c_rate}', f'c_opt : {c_opt}', f'linear_scaling : {linear_scaling}', f'function_filter : {function_filter}', f'variable_filter : {variable_filter}', f'xydomain_filter : {xydomain_filter}', f'x_domain : \n{x_domain}', f'y_domain : {y_domain}', f'equal : {equal}', '--'])

        pop_list = []
        unique_pool = []
        opt_time = 0
        
        if n_full > 0:
            n_full_ind = 0
            if txlog is not None:
                txlog.print(['full_method start'])
                _t0 = time.time()
                p10 = int(n_full/10)

            for i in range(n_full*times):
                new_ind = self._full_method(depth, node_space, c_range, c_rate, linear_scaling)
                if unique_only:
                    if new_ind.expr in unique_pool:
                        continue
                    else:
                        unique_pool.append(new_ind.expr)
                        
                if c_opt is not None:
                    opt_t0 = time.time()
                    new_ind = opt_const(new_ind, y_true=y_true, node_space=node_space, method=c_opt)
                    opt_time += time.time() - opt_t0
                    if i == 0:
                        txt_list.append('Constant node optimization... (Full)')
                        
                if function_filter or variable_filter or xydomain_filter:
                    if is_complex(  new_ind, domain_node_space_, function_filter=function_filter, 
                                    variable_filter=variable_filter, xydomain_filter=xydomain_filter,
                                    fcom=fcom, x_domain=x_domain, y_domain=y_domain, equal=equal):
                        continue

                pop_list.append(new_ind)
                n_full_ind += 1

                if n_full_ind == n_full:
                    txt_list.append(f'Number of trials (Full) : {i+1}')
                    txt_list.append(f'Number of individuals generated (Full) : {n_full_ind}')
                    if txlog is not None:
                        full_gen_end = time.time()
                        txlog.print([f'Generation time (full) : {full_gen_end-base_time}\n'])
                    break
                else:
                    if txlog is not None:
                        if p10==n_full_ind:
                            _t1 = time.time()
                            txlog.print([f'generated/trials {p10}/{i+1} (time : {_t1-_t0} s)'])
                            p10 += int(n_full/10)
                            _t0 = _t1

            if n_full != n_full_ind:
                raise Exception(f'Could not generate enough population. Generated population = {n_full_ind}. Required population = {n_full}')
        else:
            full_gen_end = time.time()

        
        if n_grow > 0:
            n_grow_ind = 0
            if txlog is not None:
                _t0 = time.time()
                txlog.print(['grow_method start'])
                p10 = int(n_grow/10)

            for i in range(n_grow*times):
                new_ind = self._grow_method(depth, node_space, c_range, c_rate, grow_50_50, linear_scaling)
                if unique_only:
                    if new_ind.expr in unique_pool:
                        continue
                    else:
                        unique_pool.append(new_ind.expr)
                
                if c_opt is not None:
                    opt_t0 = time.time()
                    new_ind = opt_const(new_ind, y_true=y_true, node_space=node_space, method=c_opt)
                    opt_time += time.time() - opt_t0
                    if i == 0:
                        txt_list.append('Constant node optimization... (Grow)')
                            
                            
                if function_filter or variable_filter or xydomain_filter:
                    if is_complex(  new_ind, domain_node_space_, function_filter=function_filter, 
                                    variable_filter=variable_filter, xydomain_filter=xydomain_filter, 
                                    fcom=fcom, x_domain=x_domain, y_domain=y_domain, equal=equal):
                        continue

                pop_list.append(new_ind)
                n_grow_ind += 1

                if n_grow_ind == n_grow:
                    txt_list.extend([f'Number of trials (Grow) : {i+1}', f'Number of individuals generated (Grow) : {n_grow_ind}\n'])
                    # txt_list.append(f'Number of individuals generated (Grow) : {n_grow_ind}\n')
                    if txlog is not None:
                        grow_gen_end = time.time()
                        txlog.print([f'Generation time (grow) : {grow_gen_end-full_gen_end}'])
                    break
                else:
                    if txlog is not None:
                        if p10==n_grow_ind:
                            _t1 = time.time()
                            txlog.print([f'generated/trials {p10}/{i+1} (time : {_t1-_t0} s)'])
                            p10 += int(n_grow/10)
                            _t0 = _t1

            if n_grow != n_grow_ind:
                raise Exception(f'Could not generate enough population. Generated population = {n_grow_ind}. Required population = {n_grow}')

        
        if txlog is not None:
            end_time = time.time()
            txlog.print([f'Generation time (all) : {end_time-base_time}\nOpt time (all) : {opt_time} s'])
        return pop_list, txt_list
    
    def _generate(self, terminal_condition, depth, node_space, c_range, c_rate, linear_scaling):
        ind = tree_graph()
        node_pool = list()
        linear_node = [dict(name='Add', arity=2, value=None), dict(name='const0', arity=0, value=0), dict(name='Mul', arity=2, value=None),  dict(name='const1', arity=0, value=1)]
        n_const = 0
        max_depth  = random.randint(*depth)
        depth_pool = [0]
        while len(depth_pool) != 0:
            current_depth = depth_pool.pop()
            if terminal_condition(current_depth, max_depth, c_rate):
                if c_rate is None:
                    c_rate = 1/(len(node_space.Xinfo) + 1)
                if c_rate > random.random():
                    node_pool.append({'name': 'const{}'.format(n_const), 'arity': 0, 'value': random.uniform(*c_range)})
                    n_const += 1
                else:
                    node_pool.append(random.choice(node_space.Xinfo))
            else:
                _node = random.choice(node_space.Finfo)
                node_pool.append(_node)
                depth_pool.extend([current_depth+1]*_node['arity'])
        if linear_scaling:
            ind.extend(linear_node)
        ind.extend(node_pool)
        ind.calc_depth()
        return ind

    def _full_method(self, depth, node_space, c_range, c_rate, linear_scaling):    
        def terminal_condition(current_depth, max_depth, c_rate):
            return current_depth == max_depth
        return self._generate(terminal_condition, depth, node_space, c_range, c_rate, linear_scaling)
    
    def _grow_method(self, depth, node_space, c_range, c_rate, grow_50_50, linear_scaling):
        def terminal_condition(current_depth, max_depth, c_rate):
            if current_depth == max_depth:
                return True
            else:
                if grow_50_50:
                    rate = 0.5 > random.random()
                else:
                    if c_rate is None:
                        c_rate = 1/(len(node_space.Xinfo) + 1)
                    num_leaf = len(node_space.Xinfo)/(1-c_rate)
                    num_func = len(node_space.Finfo)
                    rate = num_leaf / (num_leaf+num_func) > random.random()
                return (current_depth >= depth[0] and rate)
        return self._generate(terminal_condition, depth, node_space, c_range, c_rate, linear_scaling)


def select_idx_random(score_list, k):
    idx = list(range(len(score_list)))
    idx_k_time = random.choices(idx,k=k)
    return idx_k_time

def select_idx_tournament(score_list, tournament_num=5, high_best=False):
    selcted_idx = []
    for i in range(len(score_list)):
        participant_idx = select_idx_random(score_list, k=tournament_num)
        participant_score = [score_list[idx] for idx in participant_idx]
        top = participant_score.index(sorted(participant_score,reverse=high_best)[0])
        selcted_idx.append(participant_idx[top])
    return selcted_idx

def select_idx_best(score_list, high_best=False, num=1, unique=False):
    if unique:
        _, u_idx = np.unique(score_list, return_index=True)# sorted
        if high_best:
            return u_idx[::-1][:num]
        else:
            return u_idx[:num]
    else:
        top_idx = score_list.index(sorted(score_list,reverse=high_best)[:num])
        return top_idx



def Genetic_manipulation(   pop, make_tree_object, max_depth=4, p_crossover = 0.7, 
                            p_mutation = 0.1, mutat_tree_depth=(0, 2), linear_scaling=False,
                            random_select=True, detailed_log = False, node_space=None,
                            y_true = None, opt_method=None, n_trial = 1,
                            function_filter=False, variable_filter=False, xydomain_filter=False, 
                            fcom = None, x_domain=None, y_domain=None, equal=(True, True), selection=False):
    '''
    random_select : In the case of True, it is randomly selected from other than the root node. For False, select based on depth.
    '''
    detailed_log_list = ['\n', '--- Genetic manipulation detail ---']
    if linear_scaling:
        max_depth += 2
    
    def bloat_control(parent_list, child_list, max_depth):
        return_list = []
        for e, c in enumerate(child_list):
            if max(c._depth) > max_depth:
                return_list.append(parent_list[e])
            else:
                return_list.append(c)

        if detailed_log:
            if len(detailed_log_list) <= 10:
                detailed_log_list.extend(['\nparent'])
                detailed_log_list.extend([ind.expr for ind in parent_list])
                detailed_log_list.extend(['child'])
                detailed_log_list.extend([ind.expr for ind in child_list])
        return return_list
    
    def get_spoint_list(ind):
        _s = 1
        if linear_scaling: # node [0:Add, 1:const0, 2:Mul, 3:const1, 4:root of tree, 5:tree]
            _s += 4
        if random_select:
            point_1 = list(range(_s, len(ind)))
        else:
            if linear_scaling: # depth [0:Add, 1:const0, 1:Mul, 2:const1, 2:root of tree, 3:tree]
                _s -= 2
            depth_1 = random.choice(list(range(_s, max(ind._depth)+1)))
            point_1 = [idx for idx, dep in enumerate(ind._depth) if dep == depth_1]
        return point_1

    def crossover(p1, p2):
        c1, c2 = copy.deepcopy(p1), copy.deepcopy(p2)
        if linear_scaling: # node [0:Add, 1:const0, 2:Mul, 3:const1, 4:root of tree, 5:tree]
            if len(p1) < 6 or len(p2) < 6:
                return p1, p2
        else:
            if len(p1) < 2 or len(p2) <2:
                return p1, p2

        point_1 = get_spoint_list(p1)
        point_2 = get_spoint_list(p2)
        slice1 = p1.get_subtree(random.choice(point_1))
        slice2 = p2.get_subtree(random.choice(point_2))

        c1[slice1], c2[slice2] = c2[slice2], c1[slice1]
        c1.calc_depth()
        c2.calc_depth()
        return bloat_control(parent_list=[p1, p2], child_list=[c1, c2], max_depth=max_depth)
        # return bloat_control(parent_list=[p1, p2], child_list=[c1, c2], max_depth=max_depth, selection=selection)

    def mutation(p1):
        c1 = copy.deepcopy(p1)
        point_1 = get_spoint_list(p1)
        slice1 = p1.get_subtree(random.choice(point_1))
        add_tree, log_txt = make_tree_object.pop_generate(  pop_size=1, depth=mutat_tree_depth, full=0, grow=1, c_range=(-1,1), 
                                                            c_rate=None, c_opt=None, grow_50_50=False, linear_scaling=linear_scaling, 
                                                            function_filter=False, variable_filter=False, xydomain_filter=False)
        c1[slice1] = add_tree[0]
        c1.calc_depth()
        return bloat_control(parent_list=[p1], child_list=[c1], max_depth=max_depth), log_txt

    
    next_pop = [copy.deepcopy(ind) for ind in pop]

    n_crossover, n_mutation, n_try_crossover, n_try_mutation = 0, 0, 0, 0
    log_list = ['-------- Genetic_manipulation --------', '---- pop_generate (mutation) ----']
    ori_expr = [i.expr for i in next_pop]
    
    cross_time = time.time()
    for i in range(1, len(pop), 2):
        if random.random() < p_crossover:
            c1_updated, c2_updated = False, False
            for _ in range(n_trial):
                if (c1_updated and c2_updated):
                    break
                else:
                    c1, c2 = crossover(p1=next_pop[i-1], p2=next_pop[i])
                    if n_trial <= 1:
                        next_pop[i-1], next_pop[i] = c1, c2
                        c1_updated, c2_updated = True, True
                    else:
                        if c1_updated==False:
                            if c1.expr not in ori_expr:
                                next_pop[i-1], c1_updated = c1, True
                                ori_expr.append(c1.expr)
                        if c2_updated==False:
                            if c2.expr not in ori_expr:
                                next_pop[i], c2_updated = c2, True
                                ori_expr.append(c2.expr)
                n_try_crossover += 1
            n_crossover += 1
    cross_time = time.time() - cross_time
    
    ori_expr = [i.expr for i in next_pop]
    
    mutat_time = time.time()
    for i in range(len(pop)):
        if random.random() < p_mutation:
            for _ in range(n_trial):
                _ind, _log_list = mutation(p1=next_pop[i])
                if n_trial <= 1:
                    next_pop[i] = _ind[0]
                else:
                    if _ind[0].expr not in ori_expr:
                        next_pop[i] = _ind[0]
                        ori_expr.append(_ind[0].expr)
                        n_try_mutation +=1
                        break
                n_try_mutation +=1

            if n_mutation==0:
                log_list.extend(_log_list)
            n_mutation += 1
    mutat_time = time.time() - mutat_time

    opt_time = time.time()
    if opt_method is not None:
        opt_count = sum([1 for ind in next_pop if len(ind.const_list)!=0])
        opted = sum([ind.optimized for ind in next_pop])
        next_pop = [opt_const(ind, y_true=y_true, node_space=node_space, method=opt_method) for ind in next_pop]
        opted2 = sum([ind.optimized for ind in next_pop])
    else:
        opt_count = 0
        opted = 0
        opted2 = 0
    opt_time = time.time() - opt_time
    
    
    filter_time = time.time()
    for ind in next_pop:
        if is_complex(ind, node_space, function_filter=function_filter, variable_filter=variable_filter, xydomain_filter=xydomain_filter, 
                        fcom = fcom, x_domain=x_domain, y_domain=y_domain, equal=equal):
            ind.score = np.inf
    filter_time = time.time() - filter_time
            
    detailed_log_list.extend([f'\nnum Ind with const : {opt_count}\noptimized ind(Before opt) : {opted}\noptimized ind(After opt) : {opted2}',
                                f'\ntime(opt const) : {opt_time} s\ntime(crossover) : {cross_time} s\ntime(mutation) : {mutat_time} s\ntime(Apply filter) : {filter_time}',
                                f'Number of crossover : {n_crossover}', f'Number of mutation : {n_mutation}', f'Number of crossover try : {n_try_crossover}', f'Number of mutation try : {n_try_mutation}',
                                '...\n'])
    return next_pop, log_list, detailed_log_list


def opt_const(ind, y_true, node_space, method='nls'):
    if ind.const_list is not None:
        if len(ind.const_list) == 0:
            pass
        elif ind.optimized:
            pass
        else:
            if method==None:
                pass
            elif method == 'nls':
                try:
                    le_lsq=least_squares(_nls_res, x0 = ind.const_list, args=(ind, node_space, y_true), gtol=1e-5, method='lm')
                    ind.update_value(le_lsq.x)
                except:
                    pass
            else:
                minimizer_kwargs = {"args":(ind, node_space, y_true), "method":"L-BFGS-B"}
                le_lsq=basinhopping(_nls_res_b, x0 = ind.const_list, minimizer_kwargs=minimizer_kwargs, niter=10)
                ind.update_value(le_lsq.x)
    else:
        raise Exception(f'It seems that the depth of the tree has not been calculated.')
        
    return ind

def _nls_res(const, ind, node_space, y_true):
    y_pred = ind.pred(node_space.update(X=None, const_list=const))
    res = y_true - y_pred
    return res

def _nls_res_b(const, ind, node_space, y_true):
    y_pred = ind.pred(node_space.update(X=None, const_list=const))
    res = y_true - y_pred
    return sum(res**2)

def is_complex(ind, node_space, function_filter=True, variable_filter=True, xydomain_filter=True,
                fcom = [ ['sqrt'], ['Square', 'Cube'], ['log'], ['exp', 'Pow', 'Ipow'] ],
                x_domain=None, y_domain=None, equal=(True, True)):
    # 'Add', 'Sub', 'Mul', 'Div', 'Square', 'Cube', 'sqrt', 'log', 'exp', 'Pow', 'Ipow'
    com_flat = list(itertools.chain.from_iterable(fcom))
    func_pool, vals_pool = list(), list() 
    stak, old_d = 0, -1

    for e, now_d in enumerate(ind._depth):
        n = ind[e]['name']
        if old_d >= now_d:
            for _ in range(abs(old_d - now_d)+1):
                func_pool.pop()
        if function_filter:
            if  n in com_flat:
                com_bool = [n in c for c in fcom]
                consider_pairs_n = [c for b, c in zip(com_bool, fcom) if b][0]
                if sum([(n in func_pool) for n in consider_pairs_n])!=0:
                    return True
        if variable_filter: 
            if n not in node_space.Fspace:
                if n in vals_pool:
                    return True
                vals_pool.append(n)
        func_pool.append(n)
        old_d = now_d
    
    if xydomain_filter:
        if (x_domain is None or y_domain is None):
            raise NameError(f'When xydomain_filter = True, x_domain needs a dataframe and y_domain needs a tuple of minimum and maximum values.\n x_domain = {type(x_domain)}, y_domain = {y_domain}')
        y_pred = ind.pred(node_space.update(X=x_domain, const_list=ind.const_list))
        if sum(equal)==2:
            if ~np.all((y_domain[0] <= y_pred) & (y_pred <= y_domain[1])):
                return True
        elif equal[0]:
            if ~np.all((y_domain[0] <= y_pred) & (y_pred < y_domain[1])):
                return True
        elif equal[1]:
            if ~np.all((y_domain[0] < y_pred) & (y_pred <= y_domain[1])):
                return True
        else:
            if ~np.all((y_domain[0] < y_pred) & (y_pred < y_domain[1])):
                return True
    if sum([(node['value'] is None) for node in ind if node['arity']==0]) == 0:
        return True
    return False


def pop_analysis(pop, node_space):
    def count_const(ind):
        return len([1 for n in ind if 'const' in n['name']])
    
    n_pop = len(pop)

    result_dict = dict()
    for d in node_space.Finfo:
        result_dict[d['name']] = len([1 for ind in pop if d in ind])/n_pop
    for d in node_space.Xinfo:
        result_dict[d['name']] = len([1 for ind in pop if d in ind])/n_pop
    c_list = [count_const(ind) for ind in pop]    
    result_dict['const'] = sum([bool(i) for i in c_list])/n_pop
    result_dict['const/expr'] = sum(c_list)/n_pop
    return result_dict