import numpy as np
import pandas as pd
from sklearn.metrics import (max_error, mean_absolute_error,
                             mean_squared_error, r2_score)

class table_log():
    def __init__(self):
        self.log_list = []
        self.len_list = []
    def save_log(self, gen, score_list, score_name_list, print_out=True):
        log = {name:score for name, score in zip(score_name_list, score_list)}
        self.log_list.append(log)
        if print_out:
            if gen == 0:
                txt = '  Gen  '
                for name in score_name_list:
                    txt += '! -- {} -- '.format(name)
                    self.len_list.append(len(name)+8)
                self.title_txt = txt
                print(txt)
                txt = '-------'
                for l in self.len_list:
                    txt += '!' + '-'*l
                print(txt)
            txt = '{:>5}  '
            for l, s in zip(self.len_list, score_list):
                if abs(s) < 10:
                    txt += '!{:>' + '{}'.format(l-3) + '.2f}   '
                else:
                    txt += '!{:>' + '{}'.format(l-3) + '.2e}   '
            print(txt.format(gen, *score_list))
        return self.log_list
        
    def end_log(self, score_list, score_name_list, print_out=True):
        txt = '-------'
        for l in self.len_list:
            txt += '!' + '-'*l
        if print_out:    
            print(txt)
            print(self.title_txt)
        log = {name:score for name, score in zip(score_name_list, score_list)}
        self.log_list.append(log)
        txt = '-------'
        for l in self.len_list:
            txt += '!' + '-'*l
        if print_out:
            print(txt)
        txt = '{:>5}  '
        for l, s in zip(self.len_list, score_list):
            if abs(s) < 10:
                txt += '!{:>' + '{}'.format(l-3) + '.2f}   '
            else:
                txt += '!{:>' + '{}'.format(l-2) + '.2e}   '
        if print_out:
            print(txt.format('final', *score_list))
        return self.log_list


class txt_log():
    def __init__(self, file_name, save_path):
        self.file_name = file_name
        self.save_path = save_path
        self.first = True
    def print(self, txt):
        if self.first:
            mode = 'w'
            self.first = False
        else:
            mode = 'a'

        f = open(f'{self.save_path}/{self.file_name}.txt', mode)
        if isinstance(txt, list):
            pass
        elif isinstance(txt, str):
            pool = list()
            pool.append(txt)
            txt = pool
        else:
            raise Exception('input must be list or str.')
        for t in txt:
            f.write(str(t) + '\n')
        f.close()
        



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
    
    
def output_pred_val(true, pred, idx=None, e_sort=True, save_name='./004_y_pred_tr.tsv'):
    if idx is None:
        y_tr_df = dict(true=np.array(true).reshape(-1), pred=np.array(pred).reshape(-1))
    else:
        y_tr_df = dict(true=np.array(true).reshape(-1), pred=np.array(pred).reshape(-1), name=np.array(idx).reshape(-1))
    df = pd.DataFrame.from_dict(y_tr_df)
    if e_sort:
        df['error'] = df['true'] - df['pred']
        df['error'] = df['error'].abs()
        df = df.sort_values('error', ascending=False)
    df.to_csv(save_name, sep='\t')