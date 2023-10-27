import copy
import inspect
import itertools
import random
import time
import warnings
from collections import deque

import numpy as np
import pandas as pd
import sympy
from scipy import stats
from scipy.optimize import basinhopping, least_squares, leastsq, minimize

warnings.simplefilter('ignore')
from functools import partial
from multiprocessing import Pool

from joblib import Parallel, delayed
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


class NumpyBasedFunction():
    @staticmethod
    def add(x1, x2):
        return np.add(x1, x2)

    @staticmethod
    def sub(x1, x2):
        return np.subtract(x1, x2)

    @staticmethod
    def mul(x1, x2):
        return np.multiply(x1, x2)

    @staticmethod
    def div(x1, x2):
        return np.divide(x1, x2)

    @staticmethod
    def ln(x):
        return np.log(x)

    @staticmethod
    def sqrt(x):
        return np.sqrt(x)

    @staticmethod
    def square(x):
        return np.square(x)

    @staticmethod
    def cube(x):
        return np.multiply(np.square(x), x)

    @staticmethod
    def exp(x):
        return np.exp(x)

    # >>> protected functions (gplearn : https://github.com/trevorstephens/gplearn/blob/master/gplearn/functions.py)
    @staticmethod
    def protected_division(x1, x2):
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.where(np.abs(x2) > 0.001, np.divide(x1, x2), 1.)

    @staticmethod
    def protected_sqrt(x1):
        return np.sqrt(np.abs(x1))

    @staticmethod
    def protected_ln(x1):
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.where(np.abs(x1) > 0.001, np.log(np.abs(x1)), 0.)


class Node_space():
    def __init__(self, X, func_list=['add', 'sub', 'mul', 'div']):
        self._X_type = type(X)
        self._Fspace = dict()
        self._Xspace = dict()
        self._Cspace = dict()

        self._SXspace = dict()
        self._Finfo  = list()
        self._Xinfo  = list()
        self._symbol = ''

        self._set_func(func_list)
        self._set_X(X)
        self._symbol = ' '.join([i['name'] for i in self._Xinfo]) + ' CONST'
        self._set_SX(X)
        
    def Nspace(self):
        re_dict = dict()
        re_dict.update(**self._Fspace, **self._SXspace, **self._Cspace)
        return re_dict
    
    def update(self, X=None, const_list=None):
        self._set_X(X)
        self._set_C(const_list)
        return self.Nspace()
    
    def add(self, a=None, b=None): return sympy.Add(a, b)

    def sub(self, a=None, b=None): return sympy.Add(a, -b)

    def mul(self, a=None, b=None): return sympy.Mul(a, b)

    def div(self, a=None, b=None): return sympy.Mul(a, 1/b)

    def exp(self, a=None): return sympy.exp(a)

    def ln(self, a=None): return sympy.log(a)

    def square(self, a=None): return sympy.Mul(a, a)

    def cube(self, a=None): return sympy.Mul(sympy.Mul(a, a), a)
    
    def sqrt(self, a=None): return sympy.sqrt(a)
    
    def protected_division(self, a=None, b=None): return sympy.Mul(a, 1/b)
    
    def protected_sqrt(self, a=None): return sympy.sqrt(a)
    
    def protected_ln(self, a=None): return sympy.log(a)
        
    
    def _set_X(self, X=None):
        if X is not None:
            self._Xinfo = list()
            if isinstance(X, self._X_type):
                if isinstance(X, type(pd.DataFrame())):
                    self._Xspace.update({'{}'.format(x_name):np.array(X[x_name], dtype='float64') for x_name in X})
                    self._Xinfo.extend([dict(name='{}'.format(x_name), arity=0, value=None) for x_name in X])
                else:
                    npX = np.array(X)
                    nraw, ncol = npX.shape
                    self._Xspace.update({f'X{nc}':npX[:, nc] for nc in range(ncol)})
                    self._Xinfo.extend([dict(name=f'X{nc}', arity=0, value=None) for nc in range(ncol)])
            else:
                raise TypeError('The type of X does not match.')
    
    def _set_C(self, const_list=None):
        if const_list is not None:
            self._Cspace = dict()
            for e, c in enumerate(const_list):
                self._Cspace['const{}'.format(e)] = c
    
    def _set_SX(self, X):
        if isinstance(X, type(pd.DataFrame())):
            sympy.var(self._symbol)
            self._SXspace.update({x_name:eval(x_name) for x_name in X.columns})
        else:
            npX = np.array(X)
            nraw, ncol = npX.shape
            sympy.var(self._symbol)
            self._SXspace.update({'X{}'.format(nc):'X{}'.format(nc) for nc in range(ncol)})
        self._SXspace.update({'CONST':eval('CONST')})
            
    def _set_func(self, func_list=None):
        if func_list is not None:
            self._Finfo = list()
            operator_name = [m[0] for m in inspect.getmembers(self, inspect.ismethod)]
            for f_name in func_list:
                if f_name not in operator_name:
                    raise NameError('{} is not registered with the operator.'.format(f_name))
            for f in func_list:
                self._Fspace.update({f:eval('self.{}'.format(f))})
                self._Finfo.append(dict(name=f, arity=len(eval('self.{}'.format(f)).__defaults__), value=None))
    
    @property
    def Fspace(self): return self._Fspace

    @property
    def Xspace(self): return self._Xspace
    
    @property
    def Cspace(self): return self._Cspace
    
    @property
    def Finfo(self): return self._Finfo  # >>> [{'name': 'Add', 'arity': 2, 'value': None}, ..., {'name': 'Div', 'arity': 2, 'value': None}]
    
    @property
    def Xinfo(self): return self._Xinfo  # >>> [{'name': 'X0', 'arity': 0, 'value': None}, ..., {'name': 'X34', 'arity': 0, 'value': None}]
    
    @property
    def symbol(self): return self._symbol



# def add(x1, x2):
#     return np.add(x1, x2)

# def sub(x1, x2):
#     return np.subtract(x1, x2)

# def mul(x1, x2):
#     return np.multiply(x1, x2)

# def div(x1, x2):
#     return np.divide(x1, x2)

# def ln(x):
#     return np.log(x)

# def sqrt(x):
#     return np.sqrt(x)

# def square(x):
#     return np.square(x)

# def cube(x):
#     return np.multiply(np.square(x), x)

# def exp(x):
#     return np.exp(x)

# def Pow(self, a=None, b=None):
#     return np.power(a, b) if self._use_sympy else sympy.Pow(a,b)

# def Ipow(self, a=None, b=None):
#     return np.power(a, np.divide(1, b)) if self._use_sympy else sympy.Pow(a,1/b)

# # >>> protected functions (gplearn : https://github.com/trevorstephens/gplearn/blob/master/gplearn/functions.py)
# def protected_division(x1, x2):
#     with np.errstate(divide='ignore', invalid='ignore'):
#         return np.where(np.abs(x2) > 0.001, np.divide(x1, x2), 1.)


# def protected_sqrt(x1):
#     return np.sqrt(np.abs(x1))


# def protected_ln(x1):
#     with np.errstate(divide='ignore', invalid='ignore'):
#         return np.where(np.abs(x1) > 0.001, np.log(np.abs(x1)), 0.)

    # def pro_log_gplearn(self, a=None):
    #     with np.errstate(divide='ignore', invalid='ignore'):
    #         return np.where(np.abs(a) > 0.001, np.log(np.abs(a)), 0.) if self._use_sympy else sympy.log(a)

    # def pro_sqrt_gplearn(self, a=None):
    #     return np.sqrt(np.abs(a)) if self._use_sympy else sympy.sqrt(a)
    
    # def pro_div_gplearn(self, a=None, b=None):
    #     with np.errstate(divide='ignore', invalid='ignore'):
    #         return np.where(np.abs(b) > 0.001, np.divide(a, b), 1.) if self._use_sympy else sympy.Mul(a, 1/b)
    
    # def _set_func(self,func_list=None):
    #     if func_list is not None:
    #         self._Finfo = list()
    #         operator_name = [m[0] for m in inspect.getmembers(self, inspect.ismethod)]
    #         for f_name in func_list:
    #             if f_name not in operator_name:
    #                 raise NameError('{} is not registered with the operator.'.format(f_name))
    #         for f in func_list:
    #             self._Fspace.update({f:eval('self.{}'.format(f))})
    #             self._Finfo.append(dict(name=f, arity=len(eval('self.{}'.format(f)).__defaults__), value=None))