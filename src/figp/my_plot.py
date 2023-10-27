import copy
import decimal
import os

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator
from scipy.stats import rankdata
from sklearn.metrics import (max_error, mean_absolute_error,
                             mean_squared_error, mean_squared_log_error,
                             r2_score)


def gen_cmap_name(cols):#http://hydro.iis.u-tokyo.ac.jp/~akira/page/Python/contents/plot/color/colormap.html
    nmax = float(len(cols)-1)
    color_list = []
    for n, c in enumerate(cols):
        color_list.append((n/nmax, c))
    return mpl.colors.LinearSegmentedColormap.from_list('cmap', color_list)

def gen_cmap_rgb(cols, base=255):
    nmax = float(len(cols)-1)
    cdict = {'red':[], 'green':[], 'blue':[]}
    for n, c in enumerate(cols):
        loc = n/nmax
        cdict['red'  ].append((loc, c[0]/base, c[0]/base))
        cdict['green'].append((loc, c[1]/base, c[1]/base))
        cdict['blue' ].append((loc, c[2]/base, c[2]/base))
    return mpl.colors.LinearSegmentedColormap('cmap', cdict)


def narrow_cmap(cmap='viridis', start=0, stop=-20):
    x = np.linspace(0.0, 1.0, 100)
    cmap = cm.ScalarMappable(cmap=cmap)
    rgb = cmap.to_rgba(x)[np.newaxis, :, :3]
    return gen_cmap_rgb(rgb[0][slice(start,stop)], base=1.)


def lb_dr():
    # return gen_cmap_rgb([[200, 0, 0], [180, 45, 255], [145, 254, 255]][::-1])
    return gen_cmap_rgb([[200, 0, 0], [155, 45, 255], [145, 250, 255]][::-1])

def dr_lb():
    # return gen_cmap_rgb([[200, 0, 0], [180, 45, 255], [145, 254, 255]][::-1])
    return gen_cmap_rgb([[145, 250, 255], [155, 45, 255], [200, 0, 0]][::-1])

def lb_mg_dr():
    return gen_cmap_rgb([[200, 0, 0], [245, 45, 255], [145, 254, 255]][::-1])

def lb_gra_dr():
    return gen_cmap_rgb([[220, 0, 0], [135, 135, 135], [145, 224, 255]][::-1])
    

def my_cmap(name='lb_dr'):
    if name=='lb_dr':
        return lb_dr()
    elif name=='dr_lb':
        return lb_dr()
    elif name=='lb_mg_dr':
        return lb_mg_dr()
    elif name=='my_cmap':
        return gen_cmap_name(['crimson', 'magenta', 'darkviolet', 'b', 'deepskyblue', 'lime', 'g'][::-1])
    elif name=='n_viridis':
        return narrow_cmap(cmap='viridis', start=0, stop=-10)
    elif name=='n_viridis_r':
        return narrow_cmap(cmap='viridis_r', start=10, stop=None)
    elif name=='n_plasma':
        return narrow_cmap(cmap='plasma', start=0, stop=-10)
    elif name=='n_plasma_r':
        return narrow_cmap(cmap='plasma_r', start=10, stop=None)
    elif name=='n_inferno':
        return narrow_cmap(cmap='inferno', start=0, stop=-10)
    elif name=='n_inferno_r':
        return narrow_cmap(cmap='inferno_r', start=10, stop=None)
    elif name=='n_magma':
        return narrow_cmap(cmap='magma', start=0, stop=-10)
    elif name=='n_magma_r':
        return narrow_cmap(cmap='magma_r', start=10, stop=None)
    elif name=='n_cividis':
        return narrow_cmap(cmap='cividis', start=0, stop=-10)
    elif name=='n_cividis_r':
        return narrow_cmap(cmap='cividis_r', start=10, stop=None)
    else:
        return None


def _set_fig(x_minor=False, y_minor=False, font_size=30):# https://qiita.com/Miyabi1456/items/ee85861ff98c8c2c9dd0
    # plt.rcParams["font.family"] = "Times New Roman"  
    plt.rcParams["xtick.direction"] = "in"         
    plt.rcParams["ytick.direction"] = "in"          
    plt.rcParams["xtick.minor.visible"] = x_minor       
    plt.rcParams["ytick.minor.visible"] = y_minor     
    plt.rcParams["xtick.major.width"] = 2            
    plt.rcParams["ytick.major.width"] = 2          
    plt.rcParams["xtick.minor.width"] = 1.0         
    plt.rcParams["ytick.minor.width"] = 1.0          
    plt.rcParams["xtick.major.size"] = 12
    plt.rcParams["ytick.major.size"] = 12
    plt.rcParams["xtick.minor.size"] = 5            
    plt.rcParams["ytick.minor.size"] = 5            
    plt.rcParams["font.size"] = font_size                 
    plt.rcParams["axes.linewidth"] = 1.5            


def _make_legend(ax, figsize=(8,8), save_name='./legend.png', dpi=600):
    fig_leg = plt.figure(figsize=figsize)
    ax_leg = fig_leg.add_subplot(111,facecolor=(1, 1, 1, 0))
    ax_leg.legend(*ax.get_legend_handles_labels(), loc=2 ,frameon=False)
    ax_leg.axis('off')
    fig_leg.savefig(f'{save_name}', dpi=dpi, transparent=True, bbox_inches="tight", pad_inches=0.05)
    plt.clf()
    plt.close()

def _make_color_bar(cmap_obj, figsize=(3, 8), save_name='c_bar.png', dpi=600):
    fig_c_bar = plt.figure(figsize=figsize)
    ax_c_bar = fig_c_bar.add_subplot(111,facecolor=(1, 1, 1, 0))
    ax_c_bar.axis('off')
    cax = plt.axes([0.1, 0.05, 0.3/figsize[0], 0.9])
    fig_c_bar.colorbar(cmap_obj,cax=cax)
    fig_c_bar.savefig(f'{save_name}', dpi=dpi, transparent=True, bbox_inches="tight", pad_inches=0.05)
    plt.clf()
    plt.close()


def is_vector(data):
    data = np.array(data)
    # print(f'a : {data.reshape(-1).shape}\nb : {data.shape}')
    return data.reshape(-1).shape[0] == max(data.shape)

def _input_analysis(data):
    if isinstance(data, list):
        if isinstance(data[0], (float, int)):
            return '[float, ..., float]'
        else:
            if isinstance(data[0], (list, np.ndarray, pd.DataFrame, pd.Series)):
                _check_vec = [is_vector(d) for d in data]
                only_vec = np.prod(_check_vec)
                only_mat = np.sum(_check_vec)
                if only_vec==1:
                    if len(data) == 1:
                        return '[np.ndarray(vector)]'
                    else:
                        return '[np.ndarray(vector), ..., np.ndarray(vector)]'
                elif only_mat==0:
                    if len(data) == 1:
                        return '[np.ndarray(matrix)]'
                    else:
                        return '[np.ndarray(matrix), ..., np.ndarray(matrix)]'
                else:
                    raise Exception('The inputs must be aligned as vectors or matrices.')
            else:
                raise Exception('Inside list must be list, np.ndarray, pd.DataFrame, pd.Series.')

    elif isinstance(data, (np.ndarray, pd.DataFrame, pd.Series)):
        if is_vector(data):
            return 'np.ndarray(vector)'
        else:
            return 'np.ndarray(matrix)'
    else:
        raise Exception('input must be list or np.ndarray or pd.DataFrame')

def _broadcast(data, data_direction='column'):
    _dtype = _input_analysis(data)
    if _dtype == '[float, ..., float]':
        return [np.array(data)], 0
    if _dtype == '[np.ndarray(vector)]':
        return data, 0
    if _dtype == '[np.ndarray(matrix)]':
        data_list = []
        if data_direction == 'index':
            for d in data:
                d = np.array(d)
                data_list.extend([d[_d, :].reshape(-1) for _d in range(d.shape[0])])
        else:
            for d in data:
                d = np.array(d)
                data_list.extend([d[:, _d].reshape(-1) for _d in range(d.shape[1])])
        return data_list, 1
    if _dtype == '[np.ndarray(vector), ..., np.ndarray(vector)]':
        return data, 1
    if _dtype == '[np.ndarray(matrix), ..., np.ndarray(matrix)]':
        data_list = []
        if data_direction == 'index':
            for d in data:
                d = np.array(d)
                data_list.extend([d[_d, :].reshape(-1) for _d in range(d.shape[0])])
        else:
            for d in data:
                d = np.array(d)
                data_list.extend([d[:, _d].reshape(-1) for _d in range(d.shape[1])])
        return data_list, 1
    if _dtype == 'np.ndarray(vector)':
        return [np.array(data).reshape(-1)], 0
    if _dtype == 'np.ndarray(matrix)':
        data_list = []
        if data_direction == 'index':
            d = np.array(data)
            data_list.extend([d[_d, :].reshape(-1) for _d in range(d.shape[0])])
        else:
            d = np.array(data)
            data_list.extend([d[:, _d].reshape(-1) for _d in range(d.shape[1])])
        return data_list, 1
    raise Exception('Some data could not be broadcast.')

def _R24plot(true, pred, metric='R2'):
    def Significant_figures(num):
        if -10 < num <10:
            return f'{num:.3f}'
        else:
            return f'{num:.2e}'
    
    
    _bool = sum(np.isnan(true)) + sum(np.isnan(pred))
    if _bool:
        if metric == 'R2':
            return '\n$R^2$ : -inf'
        if metric == 'MAE':
            return '\nMAE : inf'
        if metric == 'RMSE':
            return '\nRMSE : inf'
    else:
        if metric == 'R2':
            score = r2_score(y_true=true, y_pred=pred)
            return f'\n$R^2$ : {Significant_figures(score)}'

        if metric == 'MAE':
            score = mean_absolute_error(y_true=true, y_pred=pred)
            return f'\nMAE : {Significant_figures(score)}'

        if metric == 'RMSE':
            score = mean_squared_error(y_true=true, y_pred=pred, squared=False)
            return f'\nRMSE : {Significant_figures(score)}'


def scatter_plot(x_data, y_data, c_data=['b', 'r', 'g'], label_data=None, xy_labels=['x', 'y'], title=None, figsize=(8,8), alpha=1, marksize=18, invert_xaxis=False,
                cmap='jet', cmin=None, cmax=None, diagonal=True, eq_aspect=False, legend=True, color_bar=True, c_bar_title='', save_name='./99_', font_size=30, dpi=600, 
                save_plain=False, return_obj=False):
    
    
    
    _set_fig(font_size=font_size)
    if label_data is None:
        label_data = [i for i in range(10000)]
    else:
        if isinstance(label_data, list):
            label_data.extend([i for i in range(10000)])
        else:
            label_data = [label_data]
            label_data.extend([i for i in range(10000)])

    _cmap = my_cmap(cmap)
    if _cmap is not None:
        cmap = _cmap
    cmap = cm.ScalarMappable(cmap=cmap)
    if isinstance(c_data, list):
        pass
    else:
        c_data = [c_data]

    if isinstance(c_data[0], str):
        color_bar = False
        c_list = c_data * 100
        c_mode = 1
    else:
        _c_list, c_mode = _broadcast(c_data, data_direction='column')
        c_list = []
        cminmax = [c for c in [cmin, cmax] if c is not None]
        if len(cminmax) == 0:
            for i in range(len(_c_list)):
                c_list.append(cmap.to_rgba(np.array(_c_list[i*c_mode])))
        else:
            for i in range(len(_c_list)):
                c_list.append(cmap.to_rgba(np.append(_c_list[i*c_mode], cminmax))[:len(cminmax)*-1])
                
    # fig = plt.figure(figsize=figsize)
    # spec = gridspec.GridSpec(ncols=2, nrows=1,
    #                          width_ratios=[1, 100])
    # ax0 = fig.add_subplot(spec[0])
    # ax0.axis("off")
    # ax = fig.add_subplot(spec[1])
    
    # fig = plt.figure(figsize=figsize)
    # spec = gridspec.GridSpec(ncols=3, nrows=3, width_ratios=[10, 100, 10], height_ratios=[10, 100, 10])
    # for _i in [0, 1, 2, 3, 5, 6, 7, 8]:
    #     exec(f'ax{_i} = fig.add_subplot(spec[_i])')
    #     eval(f'ax{_i}.axis("off")')
    #     # ax0 = fig.add_subplot(spec[_i])
    #     # ax0.axis("off")
    # ax = fig.add_subplot(spec[4])

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    
    ax.xaxis.set_major_locator(MaxNLocator(5)) 
    ax.yaxis.set_major_locator(MaxNLocator(5)) 
    ax.tick_params(axis='x', pad=float(font_size)/2)

    x_list, x_mode = _broadcast(x_data, data_direction='column')
    y_list, y_mode = _broadcast(y_data, data_direction='column')
    for i in range(max(len(x_list),len(y_list))):
        _x = x_list[i*x_mode]
        _y = y_list[i*y_mode]
        _x = np.where(np.isinf(_x), np.nan, _x)
        _y = np.where(np.isinf(_y), np.nan, _y)
        ax.scatter(_x, _y, label=str(label_data[i])+_R24plot(_x, _y), c=c_list[i*c_mode], zorder=len(_x)*-1, alpha=alpha, s=marksize)

    if diagonal:
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        ax.set_xlim(min(xmin, ymin), max(xmax, ymax))
        ax.set_ylim(min(xmin, ymin), max(xmax, ymax))
        ax.set_aspect('equal')
        plt.plot([min(xmin, ymin), max(xmax, ymax)],[min(xmin, ymin), max(xmax, ymax)], 'k-', linewidth = 1.5, zorder=0)
    elif eq_aspect:
        ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
        
    if title is not None:
        ax.set_title(title, pad=20)

    plt.xlabel(xy_labels[0])
    plt.ylabel(xy_labels[1])
    if invert_xaxis:
        ax.invert_xaxis()
    
    add = '_p'
    if return_obj:
        return fig, ax
    
    else:
        if save_plain:
            fig.savefig(f'{save_name}{add}.png', dpi=dpi, transparent=True, bbox_inches="tight", pad_inches=0.1)

        
        if legend:
            ax.legend(frameon=False, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
            add += 'l'
        if color_bar:
            fig.colorbar(cmap, label=c_bar_title)
            add += 'c'
            
        fig.savefig(f'{save_name}{add}.png', dpi=dpi ,transparent=True, bbox_inches="tight", pad_inches=0.1)
        plt.clf()
        plt.close()



def line_plot(x_data, y_data, c_data=['b', 'r', 'g'], label_data=None, xy_labels=['x', 'y'], title=None, figsize=(16,8), cmap='jet', make_legend=True, color_bar=True, c_bar_title='', cmin=None, cmax=None, save_name='./99_', data_direction='column', 
    invert_xaxis=False, linewidth = 0.5, font_size=30, facecolor=None, vspan_data=None, dpi=600, line_style=['-'], 
    save_plain=False):

    _set_fig(font_size=font_size)
    if label_data is None:
        label_data = [i for i in range(10000)]
    else:
        if isinstance(label_data, list):
            label_data.extend([i for i in range(10000)])
        else:
            label_data = [label_data]
            label_data.extend([i for i in range(10000)])
            
            
    _cmap = my_cmap(cmap)
    if _cmap is not None:
        cmap = _cmap
    cmap = cm.ScalarMappable(cmap=cmap)
    
    if isinstance(c_data, list):
        pass
    else:
        c_data = [c_data]

    if isinstance(c_data[0], str):
        color_bar = False
        c_list = c_data * 100
        c_mode = 1
    else:
        cminmax = [c for c in [cmin, cmax] if c is not None]
        
        c_data = np.array(c_data).reshape(-1)
        c_list = cmap.to_rgba(np.append(c_data, cminmax))
        c_mode = 1
        # make_legend = False
        
    fig = plt.figure(figsize=figsize)
    spec = gridspec.GridSpec(ncols=2, nrows=1,
                             width_ratios=[1, 100])
    ax0 = fig.add_subplot(spec[0])
    ax0.axis("off")
    ax = fig.add_subplot(spec[1])

    ax.xaxis.set_major_locator(MaxNLocator(5)) 
    ax.yaxis.set_major_locator(MaxNLocator(5)) 

    x_list, x_mode = _broadcast(x_data, data_direction=data_direction)
    y_list, y_mode = _broadcast(y_data, data_direction=data_direction)

    n_col = 0
    
    for i in range(max(len(x_list),len(y_list))):
        _x = x_list[i*x_mode]
        _y = y_list[i*y_mode]
        _x = np.where(np.isinf(_x), np.nan, _x)
        _y = np.where(np.isinf(_y), np.nan, _y)
        ax.plot(_x, _y, label=str(label_data[i]), c=c_list[i*c_mode], linewidth=linewidth, linestyle = line_style[i%len(line_style)])
        n_col +=1
    if vspan_data is not None:
        if isinstance(vspan_data, list):
            pass
        else:
            vspan_data = [vspan_data]
        for vspan in vspan_data:
            ax.axvspan(xmin=vspan['xmin'],xmax=vspan['xmax'], fc='r', alpha=0.1, lw=0.)

    plt.xlabel(xy_labels[0])
    plt.ylabel(xy_labels[1])
    if invert_xaxis:
        ax.invert_xaxis()
        
    if title is not None:
        ax.set_title(title)
    
    if save_plain:
        if facecolor is None:
            fig.savefig(f'{save_name}_p.png', dpi=dpi, transparent=True, bbox_inches="tight", pad_inches=0.05)
        else:
            fig.savefig(f'{save_name}_p.png', dpi=dpi, transparent=True, facecolor=facecolor, bbox_inches="tight", pad_inches=0.05)

    if make_legend or color_bar:
        _sname = 'p'
        if make_legend:
            _sname += 'l'
            # ax.legend(frameon=False, loc = 'best')
            ax.legend(frameon=False, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., ncol=(n_col//11)+1)
        if color_bar:
            _sname += 'c'
            fig.colorbar(cmap, label=c_bar_title)
        fig.savefig(f'{save_name}_{_sname}.png', dpi=dpi ,transparent=True, bbox_inches="tight", pad_inches=0.05)
    plt.clf()
    plt.close()

    # if make_legend:
    #     _make_legend(ax, save_name=f'{save_name}_l.png')
    # if color_bar:
    #     _make_color_bar(cmap, save_name=f'{save_name}_c.png')
        

def hist_plot(data, x_range=None, c_data=['b', 'r', 'g'], label_data=None, xy_labels=['x', 'y'], n_scales = (5, 5), title=None, figsize=(16,8), cmap='my_cmap', make_legend=True, save_name='./99_', data_direction='column', 
    invert_xaxis=False, font_size=30, facecolor=None, dpi=600):
    b = (n_scales[0]+1)*4
    
    _set_fig(font_size=font_size)
    if label_data is None:
        label_data = [i for i in range(10000)]
    else:
        if isinstance(label_data, list):
            label_data.extend([i for i in range(10000)])
        else:
            label_data = [label_data]
            label_data.extend([i for i in range(10000)])
            
    _cmap = my_cmap(cmap)
    if _cmap is not None:
        cmap = _cmap
    cmap = cm.ScalarMappable(cmap=cmap)
    
    if isinstance(c_data, list):
        pass
    else:
        c_data = [c_data]

    if isinstance(c_data[0], str):
        c_list = c_data * 100
        c_mode = 1
    else:
        c_data = np.array(c_data).reshape(-1)
        c_list = cmap.to_rgba(c_data)
        c_mode = 1
        
    fig = plt.figure(figsize=figsize)
    spec = gridspec.GridSpec(ncols=2, nrows=1,
                             width_ratios=[1, 100])
    ax0 = fig.add_subplot(spec[0])
    ax0.axis("off")
    ax = fig.add_subplot(spec[1])

    ax.xaxis.set_major_locator(MaxNLocator(n_scales[0])) 
    ax.yaxis.set_major_locator(MaxNLocator(n_scales[1])) 

    y_list, y_mode = _broadcast(data, data_direction=data_direction)

    n_col = 0
    for i in range(len(y_list)):
        _y = y_list[i*y_mode]
        _y = np.where(np.isinf(_y), np.nan, _y)
        ax.hist(_y, bins=b, alpha=0.5, range=x_range, color=c_list[i*c_mode], label=str(label_data[i]))
        n_col +=1

    plt.xlabel(xy_labels[0])
    plt.ylabel(xy_labels[1])
    if invert_xaxis:
        ax.invert_xaxis()
        
    if title is not None:
        ax.set_title(title)
    
    if facecolor is None:
        fig.savefig(f'{save_name}_p.png', dpi=dpi, transparent=True, bbox_inches="tight", pad_inches=0.05)
    else:
        fig.savefig(f'{save_name}_p.png', dpi=dpi, transparent=True, facecolor=facecolor, bbox_inches="tight", pad_inches=0.05)

    if make_legend:
        _sname = 'pl'
        ax.legend(frameon=False, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., ncol=(n_col//11)+1)
        fig.savefig(f'{save_name}_{_sname}.png', dpi=dpi ,transparent=True, bbox_inches="tight", pad_inches=0.05)
    plt.clf()
    plt.close()

if __name__ == '__main__':

    if 0:
        np.random.seed(seed=0)
        num = 100
        x = [np.random.rand(num), np.random.rand(num, 2), pd.DataFrame(np.random.rand(num, 2)), np.random.rand(num)]
        y = [np.random.rand(num), np.random.rand(num, 2), np.random.rand(num, 2),               [np.random.rand(num), np.random.rand(num), np.random.rand(num)]]
        for i in range(len(x)):
            scatter_plot(x_data=x[i], y_data=y[i], c_data=y[i], label_data=['train', 'test'], figsize=(8,8), font_size=30, cmap='jet', diagonal=True, save_name=f'./{i}_')
    # if 1:
    #     line_plot(x_data=Xcol, y_data=X, c_data=list(y), label_data=['train', 'test'], figsize=(16,8), font_size=30, cmap='bwr', # 'viridis', 'jet'
    #         diagonal=True, save_name=f'./098_', data_direction='index',
    #         invert_xaxis=True, linewidth=0.1, facecolor=None)
    #     line_plot(x_data=Xcol, y_data=X, c_data=list(y), label_data=['train', 'test'], figsize=(16,8), font_size=30, cmap='bwr', # 'viridis', 'jet'
    #         diagonal=True, save_name=f'./099_', data_direction='index',
    #         invert_xaxis=True, linewidth=0.1, facecolor='gainsboro')# silver, gainsboro
