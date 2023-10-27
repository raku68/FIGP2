import sys
import time
import configparser
import pandas as pd
import numpy as np
from sklearn.metrics import (max_error, mean_absolute_error,mean_squared_error, r2_score)
from sklearn.model_selection import train_test_split
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
# from PIL import Image
from IPython.display import Image
import sys
import matplotlib
matplotlib.use('Agg') 

"""for Symbilic_Reg_ST"""
import operator
import random
import copy
from deap import algorithms, base, creator, gp, tools
from figp import Symbolic_Reg
#from figp.deap_based_FGP_NLS import surveyed_individuals, FVD_filter, isfloat
#from figp.deap_based_FGP_algorithms import FGP_NLS_algorithm
#from figp.deap_based_func_node import Node_space, NumpyBasedFunction


MMS_COLUMNS = ['chembl-id', 'pot.(log,Ki)', 'pot.(nMol,Ki)', 'aromatic_smiles', 'non_stereo_aromatic_smieles',
               'all-chembl-ids', 'no.-meas.', 'pref_name', 'accession', 'natoms',
               'core', 'sub', 'sub_carbon_replacement', 'arorings', 'a_acc',
               'a_don', 'a_heavy', 'logP(o/w)', 'RBC', 'rings',
               'TPSA', 'vdw_vol', 'Weight']
MMS_COLRENAME = {"arorings": "arings", "a_acc": "acc", "a_don": "don", "logP(o/w)": "logp", "RBC": "rbc",
                 "TPSA": "tpsa", "Weight": "mw", "pot.(log,Ki)":"pot"}
                 # RBC: Rotatable Bond Counts
MMS_FEATLIST = {'10': ["arings", "acc", "don", "a_heavy", "logp", "rbc", "rings", "tpsa", "vdw_vol", "mw"],
                '7' : ["arings", "acc", "don", "logp", "rbc", "tpsa", "mw"],
                '4' : ["logp", "rbc", "tpsa", "mw"],}
MMS_PROPERTY = "pot"

config = configparser.ConfigParser()
config.read(sys.argv[1])
print(config.sections())

file       = config["FIGP"]["INPUT_FILE"]
result_dir = config["FIGP"]["RESULT_DIR"]
filter     = config["FIGP"]["FILTER"]
nfeat      = config["FIGP"]["NFEAT"] if "NFEAT" in config["FIGP"].keys() else "7"
mseparate  = config["FIGP"]["MSEPARATE"] if "MSEPARATE" in config["FIGP"].keys() else "SIMPLE"
rtrain     = float(config["FIGP"]["RTRAIN"]) if "RTRAIN" in config["FIGP"].keys() else 0.8
d_rstate   = int(config["FIGP"]["D_RSTATE"]) if "D_RSTATE" in config["FIGP"].keys() else 0
stabilize  = int(config["FIGP"]["STABILIZE"]) if "STABILIZE" in config["FIGP"].keys() else False
s_gnoise   = bool(config["FIGP"]["G_NOISE"]) if "G_NOISE" in config["FIGP"].keys() else False
s_lmd1     = float(config["FIGP"]["S_LMD1"]) if "S_LMD1" in config["FIGP"].keys() else 1.0
s_lmd2     = float(config["FIGP"]["S_LMD2"]) if "S_LMD2" in config["FIGP"].keys() else 0.5
s_clmd1    = float(config["FIGP"]["S_CLMD1"]) if "S_CLMD1" in config["FIGP"].keys() else 1.0
s_clmd2    = float(config["FIGP"]["S_CLMD2"]) if "S_CLMD2" in config["FIGP"].keys() else 0.1

# STABILIZE=1 # 0: no stabilize, 1: variable stabilize, 2: regression coefficient stabilize, 3: variable and coefficient stabilize
# S_GNOISE=0 # 0: constant noise, 1: gaussian noise
# S_LMD1=1.0 # coefficient of stability for variable (RMSE)
# S_LMD2=0.5 # magnitude of variable noise for stability check of variables
# S_CLMD1=1.0 # coefficient of stability for coefficient (RMSE)
# S_CLMD2=0.1 # magnitude of regression coefficient noise


if filter == "FVD":
    function_filter = True
    variable_filter = True 
    xydomain_filter = True
    constonly_filter= True
    domain_filter   = False
    dfilter_aug     = False
    
elif filter == "FVD2":
    function_filter = True
    variable_filter = True 
    xydomain_filter = False
    constonly_filter= True
    domain_filter   = True
    dfilter_aug     = False
    
elif filter == "FVD3":
    function_filter = True
    variable_filter = True 
    xydomain_filter = True
    constonly_filter= True
    domain_filter   = False
    dfilter_aug     = True
    
elif filter == "FVD2&3" or filter == "FVD3&2":
    function_filter = True
    variable_filter = True 
    xydomain_filter = False
    constonly_filter= True
    domain_filter   = True
    dfilter_aug     = True
    
elif filter == "FV":
    function_filter = True
    variable_filter = True 
    xydomain_filter = False
    constonly_filter= True
    domain_filter   = False
    dfilter_aug     = False

else:
    function_filter = False
    variable_filter = False
    xydomain_filter = False
    constonly_filter= True
    domain_filter   = False
    dfilter_aug     = False


mms_featlist = MMS_FEATLIST[nfeat]                
    

print("FILE:", file)
print("OUTDIR:", result_dir)
print("NFEAT:", nfeat)
print("MSEPARETE:", mseparate)
print("RTRAIN:", rtrain)
print("D_RSTATE:", d_rstate)
print("FILTER:", filter)
print("STABILIZE:", stabilize)
print("S_LMD1:", s_lmd1)
print("S_LMD2:", s_lmd2)
print("S_CLMD1:", s_clmd1)
print("S_CLMD2:", s_clmd2)
print("mms_featlist:", mms_featlist)
    
df = pd.read_table(file, index_col=0)
df = df.rename(columns=MMS_COLRENAME)
print(df.columns)
print(file, df["core"].iloc[0])
ndata = len(df.index)

X = df.loc[:, mms_featlist]
y = df.loc[:, MMS_PROPERTY]

#print(X.describe())

X_train, X_test, y_train, y_test = None, None, None, None

if mseparate == "KMEANS":
    _scaler = StandardScaler()
    _scaler.fit(X)
    X_std = _scaler.transform(X)
    #print(pd.DataFrame(X_std, columns=mms_featlist).describe())
    _kmeans = KMeans(n_clusters=2).fit_predict(X_std)
    print("kmeans:", np.sum(_kmeans == 0), np.sum(_kmeans == 1))
    _train_gp, _test_gp = (0, 1) if np.sum(_kmeans == 1) < np.sum(_kmeans == 0) else (1, 0)
    X_train, X_test, y_train, y_test = X.loc[_kmeans==_train_gp, :], X.loc[_kmeans==_test_gp, :], y.loc[_kmeans==_train_gp], y.loc[_kmeans==_test_gp]
    ntrain = len(X_train)
    #kmeans.labels_[:20]
    #print(pd.value_counts(kmeans.labels_))
    
else:
    ntrain = int(rtrain*ndata)
    X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=ntrain, random_state=d_rstate)

    
print(f"ndata: {ndata}, ntrain: {ntrain}")
print(X_train.shape, y_train.shape, X_test.shape, y_test.shape)
# ydomain = y.min(), y.max()
_range = (y.max() - y.min())
ydomain = y.min()-0.1*_range, y.max()+0.1*_range

print("output_dir", result_dir)

res = dict()
for random_state in range(5):
    print("RANDOM STATE:", random_state)
    print(f"Symblic_Reg(ST={stabilize})")
    est = Symbolic_Reg( population_size=200,
                        generations=200,
                        tournament_size=5,
                        num_elite_select=1,
                        max_depth=4,
                        function_set=('add', 'sub', 'mul', 'div', 'sqrt', 'square', 'cube', 'ln', 'exp'),
                        metric='rmse', 
                        p_crossover=0.7, 
                        p_mutation=0.2, 
                        random_state=random_state,
                        x_domain=X,
                        y_domain=ydomain,
                        var_max_trial=5000,
                        function_filter=function_filter, 
                        variable_filter=variable_filter, 
                        xydomain_filter=xydomain_filter,
                        constonly_filter=constonly_filter,
                        domain_filter   = domain_filter,
                        dfilter_aug     = dfilter_aug,
                        domain_equal    =(True, True),
                        results_dir=f"{result_dir}/{d_rstate}_{random_state}",
                        stabilize=stabilize,
                        s_gnoise=s_gnoise,
                        s_lmd1=s_lmd1,
                        s_lmd2=s_lmd2,
                        s_clmd1=s_clmd1,
                        s_clmd2=s_clmd2)


    # traininig
    est.fit(X_train, y_train)
    score_min = float(pd.DataFrame(est.log)["score-min"].min())
    print(score_min)
    y_train_pred = est.predict(X_train)
    rmse_train = mean_squared_error(y_true=y_train, y_pred=y_train_pred, squared=False)
    r2_train = r2_score(y_true=y_train, y_pred=y_train_pred)
    res[random_state] = (score_min, rmse_train, r2_train, est)
    est.save_all()

# save the training results
res_vals = sorted(res.values(), key=lambda x: (x[0],x[1]))
print("train all results: score, rmse, r2")
print([(val[0], val[1], val[2]) for val in res_vals])
print("TRAIN RESULTS: RMSE, R2")
#print((res_vals[0][0],res_vals[0][1])) # **bug!!** fixed 2023/03/20. Before then, read results from "train all results" .
print((res_vals[0][1],res_vals[0][2]))
best_model = res_vals[0][3]
#best_model.save_all()

y_test_pred = best_model.predict(X_test)
rmse_test = mean_squared_error(y_true=y_test, y_pred=y_test_pred, squared=False)
r2_test = r2_score(y_true=y_test, y_pred=y_test_pred)

print("TEST RESULTS: RMSE, R2")
print((rmse_test, r2_test))
