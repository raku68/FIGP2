# FIGP2 for QSAR/QSPR
FIGP2 represents an enhancement in the generalizability over its predecessor, FIGP (https://github.com/takakikatsushi/FIGP). FIGP employs a symbolic regression approach, utilizing Filter-introduced Genetic Programming coupled with constant optimization via Nonlinear Least Square. The implementation of Genetic Programming is facilitated by the **DEAP** library, while the **sympy** library is employed for the analysis of mathematical expressions. Four filters can be applied as expression filters: Variable filter, Functional filter, Domain filter, and Domain filter 2.

The main target of FIGP2 is for QSAR/QSPR models. The three filters help find human-understandable mathematical expressions as a result of giving up the extremely precise expressions explaining a training data set. Details are found in our publication. 


### Reference: Shirasawa R.; Takaki K.; Miyao T. Generalizability Improvement of Symbolic Regression Models for Interpretable Quantitative Structure-activity Relationship Models, in preparation

## Gettingb Started
### Prerequisites
The following libraries are necessary on top of [Python 3.9](https://www.python.org/downloads/release/python-390/).

* [DEAP](https://github.com/DEAP/deap)  
* [scipy](https://github.com/scipy/scipy)  
* [sympy](https://github.com/sympy/sympy)  
* [matplotlib](https://github.com/matplotlib/matplotlib)
* [pandas](https://github.com/pandas-dev/pandas)
* [scikit-learn](https://github.com/scikit-learn/scikit-learn)

  
### Installation
We recomend that above packages are installed before running the installation commnad, although the installation command can automatically install theose packages.


1. Close the github repository
```
git clone https://github.com/raku68/FIGP2.git 
```

2. Create virtual environment (e.g. conda)
```
conda create -n figp2 python=3.9
conda activate figp2
```

3. Then, move to the main folder (FIGP) and install the library
```
python setup.py install
```

4. Check whether the library is installed or not.
```
$> python
Python 3.9.13 (main, Aug 25 2022, 18:29:29)
[Clang 12.0.0 ] :: Anaconda, Inc. on darwin
Type "help", "copyright", "credits" or "license" for more information.
>>> from figp import Symbolic_Reg
>>> Symbolic_Reg
<class 'figp.deap_based_FGP_NLS.Symbolic_Reg'>
```

## How to use  
### Symbolic Regression
Running symbolic regression is straightforward. All you need X and y, which are `Pandas.DataFrame` and `Pandas.Series`, respectively. 
```
from figp import Symbolic_Reg
est = Symbolic_Reg( population_size=200,
                    generations=100,
                    x_domain=X,
                    y_domain=(0, 1),
                    results_dir='./result'
                    )
est.fit(X_train, y_train)
est.save_all()
```
Input arguments (parameters) are explained in the source codes.

### Regression outcomes
The outcomes of the symbolic regression are stored in the directory specified by `results_dir` argument. In the example code above, the `result` folder is created and files are stored in the folder.
The folder contains the following files.
```
./result/000_GP_log_txt.txt',
 './result/001_FilterGPSR_params.tsv',
 './result/001_GP_log.tsv',
 './result/001_GP_log_min_pl.png',
 './result/001_GP_log_pl.png',
 './result/001_GP_log_std_pl.png',
 './result/001_GP_log_unique_rate_pl.png',
 './result/001_GP_node_analysis.tsv',
 './result/001_GP_node_analysis_func_all_pl.png',
 './result/001_GP_node_analysis_func_all_select_pl.png',
 './result/001_GP_node_analysis_func_pl.png',
 './result/001_GP_node_analysis_func_select_pl.png',
 './result/001_GP_node_analysis_select.tsv',
 './result/001_GP_node_analysis_X_pl.png',
 './result/001_GP_node_analysis_X_select_pl.png',
 './result/002_best_expr4word.txt',
 './result/002_best_expr4word_3digits.txt',
 './result/002_best_expr4word_3digits_non_e.txt',
 './result/002_best_expr4word_non_e.txt',
 './result/002_best_model.tex',
 './result/002_best_model_expr.png']
```
All the expressions searched during evolution can be found in `001_GP_log.tsv`. You can monitor tested expressions from this file.
1. Files with the prefix of **001_GP_log** store RMSE values (fitness scores) during generation.  
2. Files with the prefix of **001_GP_node** show node profiles during generation.  
3. Files with the prefix of **002_best** contain the information of best formula.
   

The best expression can be stored as an human readable expression: `002_best_model_expr.png`

![result pic](tmp/002_best_model_expr.png)

**Note!!** the variable names in `Pandas.DataFrame` containing **hyphens** and **spaces** raise Exception. Furthermore, **underscores** are recognized as subscripts when translating to expressions. Our recommendation is to use CamelStyle for variable names.

## Authors 
Raku Shirasawa: https://github.com/raku68
Takaki Katsushi: https://github.com/takakikatsushi

### Contributors to the FIGP descriptors project:
Raku Shirasawa
Takaki Katsushi
Tomoyuki Miyao

## License
See LICENSE.txt
