import os

import numpy as np
import pandas as pd

import pickle

# %matplotlib inline
from bokeh.io import show
from bokeh.layouts import column
from bokeh.models import RangeTool
from bokeh.plotting import figure, show, output_file, save
from bokeh.palettes import Magma, Inferno, Plasma, Viridis, Cividis, Colorblind, Bokeh

def read_df(file_name, compression_opts):
    try:
        df = pd.read_csv(file_name, compression=dict(method='zip', archive_name='out.csv'), index_col=0)
    except:
        df = pd.read_csv(file_name, compression=compression_opts, index_col=0)
    df.index = pd.to_datetime(df.index)
    return df

class PlotANN(object):
    """
    A class to plot datasets from ANN predictions using the ModelANN object
	Default stations are: 'Rock Slough','Jersey Point',	'Emmaton',	'Collinsville',	'Chipps',	'San Andreas Landing'
    """
    def __init__(self, dsp_home, experiments, save_dir,
				 stations = ['ROLD024', 'RSAN018', 'RSAC092', 'RSAC081', 'RSAC075', 'RSAN032','RSAN007','SLMZU011'],
				 model_file="mtl_i90_lstm14_lstm14_f_o1.h5",
                 compression_opts=None, #or dict(method='zip', archive_name='out.csv')
				 local_root_path =  './'):
                 
        self.dsp_home = dsp_home
        self.experiments = experiments
        self.save_dir = save_dir
        self.stations = stations
        self.model_file = model_file
        self.compression_opts =  compression_opts
        self.local_root_path = local_root_path

        self.col_num = len(self.stations)
        
        self.model_name = os.path.splitext(model_file)[0]

        if not os.path.exists(self.save_dir): os.makedirs(self.save_dir)

        self.exp_keys = [exp+f"_{self.model_name}" for exp in self.experiments]
        self.col_names = [e+"_EC" for e in self.experiments]
    
    def store_tuples(self,
                     file_name = "dsm2_ann_inputs_historical.csv"):
        # this function is very fast
        
        file_part = os.path.splitext(file_name)[0]

        plot_tuples = []

        if self.experiments[0] is not None:
            experiment_dir = os.path.join(self.local_root_path, "Experiments", self.experiments[0])
            target_file = os.path.join(experiment_dir, "results", "target", file_part + "_target.csv")
            target = read_df(target_file, self.compression_opts)
            plot_tuples.append((target, target.columns[self.col_num], self.experiments[0] + "_target"))

        for experiment in self.experiments:
            experiment_dir = os.path.join(self.local_root_path, "Experiments", experiment)

            model_prediction_dir = os.path.join(experiment_dir, "results", "prediction", self.model_name)
            prediction_file = os.path.join(model_prediction_dir, file_name)
            if not os.path.exists(prediction_file):
                raise ValueError(f"prediction_file {prediction_file} does not exist")
            else:
                prediction = read_df(prediction_file, self.compression_opts)
                plot_tuples.append((prediction, prediction.columns[self.col_num], experiment + "_" + self.model_name))

        self.plot_tuples = plot_tuples
			
            # Create roving RMSE windows to compare RMSE from colab versus 6years.

    def load_RMSE(self,rmse_win_len=10,rmse_calced=False):
        if rmse_calced:
            print('Reading RMSE')
            with open(f'rmse_{".".join(self.exp_keys)}.pickle', 'rb') as handle:
                self.plt_rmse_tup = pickle.load(handle)

        else:
            self.RMSE_window(rmse_win_len=rmse_win_len)

    def RMSE_window(self,
                    rmse_win_len=10):

        # this function takes up to 15 minutes
        # rmse_win_len is the window length (default 10 days) within which to create a roving RMSE value

        target_tup = self.plot_tuples[0]
        if 'target' not in target_tup[2]:
            raise(f'First tuple element needs to be "target", not {target_tup[2]}')
        
        experiment_names = [d[2] for d in self.plot_tuples]

        rmse_win_tup = {}

        for s in target_tup[0].columns:
            if s in self.stations:
                rmse_win_df = pd.DataFrame(index=target_tup[0].index[:-rmse_win_len], columns=experiment_names[1:])
                
                print(f'calculating station {s}')

                for t, tuple in enumerate(self.plot_tuples[1:]):
                    exper_tup = tuple
                    print(f'Experiment to be compared is {exper_tup[2]}')
                    print(f'Experiment {exper_tup[2]} first s of station {s}: {exper_tup[0][s][0:3].values}')

                    for i in rmse_win_df.index:
                        target_days = target_tup[0][s][i:i+pd.Timedelta(days=rmse_win_len)]
                        rmse_win_df[exper_tup[2]][i] = np.sqrt(((exper_tup[0][s][i:i+pd.Timedelta(days=10)]-target_days)**2).mean())

                rmse_win_tup[s] = rmse_win_df

        # Store the RMSE tuple 
        self.rmse_win_tup = rmse_win_tup
        
        # Export the RMSE dataframes to find correlations with other data
        out_col_names = ['target_EC']
        out_col_names.extend([exp+"_EC" for exp in experiments])
        out_col_names.extend([exp+"_RMSE" for exp in experiments])

        plt_tup_order = []
        for pt in self.plot_tuples[1:]: # don't include target
            plt_tup_order.append(pt[2])

        plt_rmse_tup = {}

        for station in self.stations:
            out_df = pd.DataFrame(index=rmse_win_tup[station].index, columns=out_col_names)

            out_df['target_EC'] = self.plot_tuples[0][0][station]

            for e, exp in enumerate(experiments):
                print(f'Checking that {exp+"_EC"} corresponds to {plt_tup_order[e]}')
                out_df[exp+"_EC"] = self.plot_tuples[e+1][0][station]

                print(f'Checking that {exp+"_RMSE"} corresponds to {self.exp_keys[e]}')
                out_df[exp+"_RMSE"] = self.rmse_win_tup[station][self.exp_keys[e]]

            plt_rmse_tup[station] = out_df

        # Store the RMSE plot tuple
        self.plt_rmse_tup = plt_rmse_tup
        with open(f'rmse_{".".join(self.exp_keys)}.pickle', 'wb') as handle:
            pickle.dump(plt_rmse_tup, handle, protocol=pickle.HIGHEST_PROTOCOL)

    def store_plot(self,
                  sta_name,
                  colors = ['#d62728', '#1f77b4', '#ff7f0e', '#2ca02c', '#9467bd', 
                            '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']): 
        # default colors same as sns.color_palette("tab10", n_colors=10).as_hex() but starts with red
        
        plt_df = self.plt_rmse_tup[station]
        
        # get the first DateTime in the dataframe index
        first = plt_df.index[0]
        # get the last DateTime in the dataframe index
        last = plt_df.index[-1]

        p = figure(height=250, width=900, tools = "xpan,wheel_zoom,box_zoom,reset,save,hover",
                x_axis_type="datetime", x_axis_location="above",
                background_fill_color="#efefef", x_range=(first, last), title=f'{sta_name} EC data')

        p.line(plt_df.index, plt_df['target_EC'], color='black', legend_label='Target')

        for i, exp in enumerate(self.experiments):
            p.line(plt_df.index, plt_df[exp+"_EC"], color=colors[i], legend_label=exp)

        r = figure(height=250, width=900, tools = "xpan,wheel_zoom,box_zoom,reset,save,hover",
                x_axis_type="datetime", x_axis_location="above", x_range=p.x_range,
                background_fill_color="#efefef", title=f'{sta_name} RMSE')

        for i, col_name in enumerate(self.experiments):
            r.line(plt_df.index, plt_df[col_name+"_RMSE"], color=colors[i], legend_label=col_name)

        select = figure(title="Drag the middle and edges of the selection box to change the range",
                        height=130, width=900, y_range=p.y_range,
                        x_axis_type="datetime", y_axis_type=None,
                        tools="", toolbar_location=None, background_fill_color="#efefef")

        range_tool = RangeTool(x_range=p.x_range)
        range_tool.overlay.fill_color = "navy"
        range_tool.overlay.fill_alpha = 0.2

        for i, col_name in enumerate(self.experiments):
            select.line(plt_df.index, plt_df[col_name+"_EC"], color=colors[i])

        select.ygrid.grid_line_color = None
        select.add_tools(range_tool)

        diff = figure(title=f"Difference (Experiment - {self.experiments[0]})",
                    height=200, width=900,
                    x_axis_type="datetime",  x_range=p.x_range,
                    tools="", toolbar_location=None, background_fill_color="#efefef")
        diff.line(plt_df.index, 0, color='black')
        
        for i, col_name in enumerate(self.experiments[1:]):
            diff.line(plt_df.index, plt_df[col_name+"_EC"] - plt_df[self.experiments[0]+"_EC"], color=colors[1+i])
        
        diff_rmse = figure(title=f"RMSE Difference (Experiment - {self.experiments[0]})",
                    height=200, width=900,
                    x_axis_type="datetime",  x_range=p.x_range,
                    tools="", toolbar_location=None, background_fill_color="#efefef")
        diff_rmse.line(plt_df.index, 0, color='black')

        for i, col_name in enumerate(self.experiments[1:]):
            diff_rmse.line(plt_df.index, plt_df[col_name+"_RMSE"] - plt_df[self.experiments[0]+"_RMSE"], color=colors[1+i])

        plot_final = column(p,r, select, diff, diff_rmse)
        if self.save_dir is not None:
            output_file(f'{self.save_dir}/{sta_name}.html', mode='inline')
            save(plot_final)
        else:
            show(plot_final)
	

if __name__ == '__main__':
	
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    dsp_home = "./"
    experiments = ["latinhypercube_v3"]
    save_dir = 'plots/lathypcub_v3'

    # initialize PlotANN object
    annplt = PlotANN(dsp_home, experiments, save_dir,model_file="mtl_i90_lstm12_lstm12_f_o1.h5")

    # create plot tuples
    annplt.store_tuples()

    # create roving RMSE tuples
    annplt.load_RMSE(rmse_calced=False)

    # create plots 
    for station in annplt.stations:
        print(f'Making plot for station {station}')
        annplt.store_plot(station)