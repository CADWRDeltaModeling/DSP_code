import os
import sys
import time

import pandas as pd
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tqdm.auto import tqdm, trange

from datetime import datetime

tqdm.pandas()

import joblib


os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
tf.config.list_physical_devices('GPU')
global local_root_path
local_root_path = "."
sys.path.append(local_root_path)
import annutils
import importlib
importlib.reload(annutils)

# Functions
def mse_loss_masked(y_true, y_pred):
    squared_diff = tf.reduce_sum(tf.math.squared_difference(y_pred[y_true > 0], y_true[y_true > 0]))
    return squared_diff / (tf.reduce_sum(tf.cast(y_true > 0, tf.float32)) + 0.01)

def predict(location, df_input, output_columns):
    model=keras.models.load_model('%s.h5'%location,custom_objects={"mse_loss_masked": mse_loss_masked})
    xscaler,yscaler=joblib.load('%s_xyscaler.dump'%location)
    return predict_with_model(model, xscaler, yscaler, df_input, output_columns)

def predict_with_model(model, xscaler, yscaler, df_input, output_columns):
    dfx = pd.DataFrame(xscaler.transform(df_input), df_input.index, columns=df_input.columns)

    yyp=model.predict(dfx, verbose=True)
    predicted_y = yscaler.inverse_transform(yyp)
    return pd.DataFrame(predicted_y, index=df_input.index, columns=output_columns)

def basic_1d(
        filters,
        stage=0,
        block=0,
        kernel_size=3,
        numerical_name=False,
        stride=None,
        force_identity_shortcut=False
):
    """
    A one-dimensional basic block.
    :param filters: the output’s feature space
    :param stage: int representing the stage of this block (starting from 0)
    :param block: int representing this block (starting from 0)
    :param kernel_size: size of the kernel
    :param numerical_name: if true, uses numbers to represent blocks instead of chars (ResNet{101, 152, 200})
    :param stride: int representing the stride used in the shortcut and the first conv layer, default derives stride from block id
    """
    parameters = {
        "kernel_initializer": "he_normal"
        }
    
    if stride is None:
        if block != 0 or stage == 0:
            stride = 1
        else:
            stride = 2

    if block > 0 and numerical_name:
        block_char = "b{}".format(block)
    else:
        block_char = chr(ord('a') + block)

    stage_char = str(stage + 2)

    def f(x):
        y = keras.layers.ZeroPadding1D(padding=1, name="padding{}{}_branch2a".format(stage_char, block_char))(x)
        y = keras.layers.Conv1D(filters, kernel_size, strides=stride, use_bias=False,
                                name="res{}{}_branch2a".format(stage_char, block_char),
                                **parameters)(y)
        y = keras.layers.BatchNormalization()(y)
        y = keras.layers.Activation("relu", name="res{}{}_branch2a_relu".format(stage_char, block_char))(y)

        y = keras.layers.ZeroPadding1D(padding=1, name="padding{}{}_branch2b".format(stage_char, block_char))(y)
        y = keras.layers.Conv1D(filters, kernel_size, use_bias=False,
                                name="res{}{}_branch2b".format(stage_char, block_char),
                                **parameters)(y)
        y = keras.layers.BatchNormalization()(y)

        if block != 0 or force_identity_shortcut:
            shortcut = x
        else:
            shortcut = keras.layers.Conv1D(filters, 1, strides=stride, use_bias=False,
                                           name="res{}{}_branch1".format(stage_char, block_char),
                                           **parameters)(x)
            shortcut = keras.layers.BatchNormalization()(shortcut)

        y = keras.layers.Add(name="res{}{}".format(stage_char, block_char))([y, shortcut])

        y = keras.layers.Activation("relu", name="res{}{}_relu".format(stage_char, block_char))(y)

        return y

    return f

def build_layer_from_string_def(s='i120', width_multiplier=1,
                                block=0,
                                force_identity_shortcut=False,
                                return_sequences_rnn=True):

    if s[0:4] == 'lstm':
        units=int(s[4:]) 
        print("building with return_sequences ",return_sequences_rnn)
        retval = layers.LSTM(units, return_sequences=return_sequences_rnn,  # was 12
                           activation='sigmoid')
        return retval

    elif s[0:3] == 'res':
        fields = s[3:].split('x')
        return basic_1d(filters=int(fields[0]),
                        stage=int(fields[3]),
                        block=block,
                        kernel_size=int(fields[1]),
                        stride=int(fields[2]),
                        force_identity_shortcut=force_identity_shortcut)
    elif s[0:3] == 'c1d':
        fields = s[3:].split('x')
        return keras.layers.Conv1D(filters=int(fields[0]), kernel_size=int(fields[1]), strides=int(fields[2]),
                                   padding='causal', activation='linear')
    elif s[0:2] == 'td':
        return keras.layers.TimeDistributed(keras.layers.Dense(int(s[2:]), activation='elu'))
    elif s[0:2] == 'dr':
        return keras.layers.Dropout(float(s[2:]))
    elif s[0] == 'f':
        return keras.layers.Flatten()
    elif s[0] == 'g':
        return keras.layers.GRU(int(s[1:]) * width_multiplier, return_sequences=True, activation='relu')
    elif s[0] == 'd':
        return keras.layers.Dense(int(s[1:]) * width_multiplier, activation='elu')
    elif s[0] == 'o':
        return keras.layers.Dense(int(s[1:]) * width_multiplier, activation='linear')
    else:
        raise Exception('Unknown layer def: %s' % s)

def build_model_string(model_type, num_neurons_multiplier, input_dropout=0, intermediate_dropout=0):
    model_type = model_type.lower()
    model_str_def = None
    if model_type == 'lstm':
        # 2. LSTM Network
        model_str_def = '%slstm%d_lstm%d_%sf_o1' % (('dr%.2f_' % input_dropout if input_dropout > 0 else ''),
                                             num_neurons_multiplier[0],num_neurons_multiplier[0],
                                             ('dr%.2f_' % intermediate_dropout if intermediate_dropout > 0 else ''),)

    return model_str_def

class ModelANN(object):
    """
    A class to create datasets to input to ANN
    """
    def __init__(self, dsp_home, 
                 compression_opts=None, #or dict(method='zip', archive_name='out.csv')
                 ndays=118, NFEATURES=8, initial_lr = 0.01,
                 stas_include=None):
                 
        self.dsp_home = dsp_home
        self.compression_opts =  compression_opts
        self.ndays = ndays
        self.NFEATURES = NFEATURES
        self.initial_lr = initial_lr

        now = datetime.now()
        root_logdir = os.path.join(os.curdir, "tf_training_logs", now.strftime("%Y%m%d-%H%M%S"))

        self.tensorboard_cb = keras.callbacks.TensorBoard(root_logdir)

        # Make a dir named Experiments
        if not os.path.exists("Experiments"):
            os.mkdir("Experiments")

        self.num_sheets = 9
        
        if stas_include is not None:
            self.observed_stations_ordered_by_median = stas_include
        else:
            self.observed_stations_ordered_by_median = ['RSMKL008', 'RSAN032', 'RSAN037', 'RSAC092', 'SLTRM004', 'ROLD024',
                                                'CHVCT000', 'RSAN018', 'CHSWP003', 'CHDMC006', 'SLDUT007', 'RSAN072',
                                                'OLD_MID', 'RSAN058', 'ROLD059', 'RSAN007', 'RSAC081', 'SLMZU025',
                                                'RSAC075', 'SLMZU011', 'SLSUS012', 'SLCBN002', 'RSAC064']

        self.output_stations = ['CHDMC006-CVP INTAKE', 'CHSWP003-CCFB_INTAKE', 'CHVCT000-VICTORIA INTAKE',
                        'OLD_MID-OLD RIVER NEAR MIDDLE RIVER', 'ROLD024-OLD RIVER AT BACON ISLAND',
                        'ROLD059-OLD RIVER AT TRACY BLVD', 'RSAC064-SACRAMENTO R AT PORT CHICAGO',
                        'RSAC075-MALLARDISLAND', 'RSAC081-COLLINSVILLE', 'RSAC092-EMMATON',
                        'RSAC101-SACRAMENTO R AT RIO VISTA', 'RSAN007-ANTIOCH', 'RSAN018-JERSEYPOINT',
                        'RSAN032-SACRAMENTO R AT SAN ANDREAS LANDING', 'RSAN037-SAN JOAQUIN R AT PRISONERS POINT',
                        'RSAN058-ROUGH AND READY ISLAND', 'RSAN072-SAN JOAQUIN R AT BRANDT BRIDGE',
                        'RSMKL008-S FORK MOKELUMNE AT TERMINOUS', 'SLCBN002-CHADBOURNE SLOUGH NR SUNRISE DUCK CLUB',
                        'SLDUT007-DUTCH SLOUGH', 'SLMZU011-MONTEZUMA SL AT BELDONS LANDING',
                        'SLMZU025-MONTEZUMA SL AT NATIONAL STEEL', 'SLSUS012-SUISUN SL NEAR VOLANTI SL',
                        'SLTRM004-THREE MILE SLOUGH NR SAN JOAQUIN R', 'SSS-STEAMBOAT SL', 'CCW-MIDDLE RIVER INTAKE',
                        'OH4-OLD R @ HWY 4', 'SLRCK005-CCWD_Rock', 'MRU-MIDDLE RIVER AT UNDINE ROAD', 'HLL-HOLLAND TRACT',
                        'BET-PIPER SLOUGH @ BETHEL TRACT', 'GES-SACRAMENTO R BELOW GEORGIANA SLOUGH',
                        'NMR: N FORK MOKELUMNE R NEAR WALNUT GROVE', 'IBS-CORDELIA SLOUGH @ IBIS CLUB',
                        'GYS-GOODYEAR SLOUGH AT MORROW ISLAND CLUB', 'BKS-SLBAR002-North Bay Aqueduct/Barker Sl']

        self.output_stations, self.name_mapping = annutils.read_output_stations(self.output_stations, self.observed_stations_ordered_by_median)
        print(self.output_stations)
        print(self.name_mapping)
    
    def build_model_from_string_def(self, strdef='i120_f_d4_d2_d1', width_multiplier=None,units=None):
        if width_multiplier is None and units is None:
            raise ValueError("width_multiplier or units must be supplied")
        print("width_multiplier")
        layer_strings = strdef.split('_')
        print ('layer_strings:%s' % layer_strings)
        inputs = keras.layers.Input(shape=[int(layer_strings[0][1:]) * self.NFEATURES])
        x = None
        prev_conv_output_num_of_channels = None
        return_sequences_rnn = True
        for block, f in enumerate(layer_strings[1:-1]):
            if x is None:
                if ('lstm' in strdef) or ('g' in strdef):
                    # these layers require 2D inputs and permutation
                    print("yo",return_sequences_rnn)
                    x = layers.Reshape((self.ndays + self.nwindows, self.NFEATURES))(inputs)
                    prev_conv_output_num_of_channels = self.NFEATURES
                    #x = layers.Permute((2, 1))(x)
                    #return_sequences_rnn = layer_strings[block + 2].startswith(('lstm', 'g', 'res', 'c1d'))
                elif ('res' in strdef) or ('cld' in strdef):
                    # these layers require 2D inputs
                    x = layers.Reshape((self.ndays + self.nwindows, self.NFEATURES))(inputs)
                    prev_conv_output_num_of_channels = self.NFEATURES
                else:
                    x = inputs

            x = build_layer_from_string_def(f, width_multiplier, block,
                                            force_identity_shortcut=(
                                                    f.startswith('res') and prev_conv_output_num_of_channels == int(
                                                f[3:].split('x')[0])),
                                            return_sequences_rnn=return_sequences_rnn)(x)
            return_sequences_rnn=False
            if f.startswith('lstm'):
                prev_conv_output_num_of_channels = int(f[4:])
            elif f.startswith('res') or f.startswith('c1d'):
                prev_conv_output_num_of_channels = int(f[3:].split('x')[0])

        out_dense_size = 23 #int(layer_strings[-1][1:]) * width_multiplier
        print(f"Output dense layer size: {out_dense_size}")
        outputs = keras.layers.Dense(out_dense_size, activation='linear')(x)
        model = keras.Model(inputs=inputs, outputs=outputs)
        model.compile(optimizer=keras.optimizers.Adam(
            learning_rate=self.initial_lr), loss="mse")
        return model

    def build_or_load_model(self, model_path, model_str_def, output_shape):
        xscaler = None
        yscaler = None
        if os.path.exists(model_path + '.h5'):
            loaded_model = annutils.load_model(model_path,
                                            custom_objects={"mse_loss_masked": mse_loss_masked})
            model = loaded_model.model
            xscaler = loaded_model.xscaler
            yscaler = loaded_model.yscaler
            print('Ignored defined model arc and loaded pre-trained model from %s.h5' % model_path)

        len_stations = output_shape[1]
        print("len_stations: ", len_stations)
        if 'lstm' in model_str_def.lower():
            model = self.build_model_from_string_def(model_str_def, width_multiplier=len_stations)

        return model, xscaler, yscaler

    def concat_dfs(self, files, windows):

        X_df= None
        Y_df= None
        for infile, window in zip(tqdm(files), windows):
            if not os.path.exists(infile):
                print("local root: ",local_root_path)
                raise ValueError(f"Data path {infile} does not exist")
            dfinps, dfouts = annutils.read_and_split(infile, self.num_sheets, self.observed_stations_ordered_by_median, vars_include=self.vars_include)
            print("\nread_split",dfinps.first_valid_index())
            
            dfinps = annutils.create_antecedent_inputs(dfinps,ndays=self.ndays,window_size=self.window_size,nwindows=self.nwindows)
            dfinps, dfouts = annutils.synchronize(dfinps, dfouts)
            
            dfinps = annutils.include(dfinps, window)
            dfouts = annutils.include(dfouts, window)
            X_df = pd.concat([X_df, dfinps], axis=0)
            Y_df = pd.concat([Y_df, dfouts], axis=0)

        return X_df, Y_df
    
    def compile_inputs(self, experiment, train_files, train_windows, test_files, test_windows, 
                       window_size=0, nwindows=0, vars_include=None):
        
        self.window_size = window_size
        self.nwindows = nwindows
        self.experiment = experiment
        self.vars_include = vars_include

        if not os.path.exists("Experiments/" + experiment):
            os.mkdir("Experiments/" + experiment)

        train_X, train_Y = self.concat_dfs(train_files, train_windows) 
        test_X, test_Y = self.concat_dfs(test_files, test_windows) 

        xx= list(train_X.columns.copy())
        yy = list(test_X.columns.copy())

        for xxx,yyy in zip(xx,yy): 
            if "dcc" in xxx or "dcc" in yyy or "128" in xxx or "128" in yyy: 
                print(xxx," ",yyy)
        nonovlp = len(set(train_X.columns) ^ set(test_X.columns))
        if nonovlp > 0:
            raise ValueError(f'Non-overlapping columns in train and test data: {nonovlp}')

        train_X.to_csv(os.path.join("Experiments", experiment, "train_X.csv"), compression=self.compression_opts,float_format="%.2f")
        train_Y.to_csv(os.path.join("Experiments", experiment, "train_Y.csv"), compression=self.compression_opts,float_format="%.2f")
        test_X.to_csv(os.path.join("Experiments", experiment, "test_X.csv"), compression=self.compression_opts,float_format="%.2f")
        test_Y.to_csv(os.path.join("Experiments", experiment, "test_Y.csv"), compression=self.compression_opts,float_format="%.2f")

        print(f"Finished compiling inputs for {experiment} experiment")

    def train_models(self, models):
        print("experiment: ", self.experiment)

        # create folders to save results
        result_folders = ['models', 'results', 'images']
        for result_folder in result_folders:
            folder_path = os.path.join(local_root_path, "Experiments", self.experiment, result_folder)
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)

        train_X = pd.read_csv(os.path.join("Experiments", self.experiment, "train_X.csv"), index_col=0, compression=self.compression_opts)
        train_Y = pd.read_csv(os.path.join("Experiments", self.experiment, "train_Y.csv"), index_col=0, compression=self.compression_opts)
        test_X = pd.read_csv(os.path.join("Experiments", self.experiment, "test_X.csv"), index_col=0, compression=self.compression_opts)
        test_Y = pd.read_csv(os.path.join("Experiments", self.experiment, "test_Y.csv"), index_col=0, compression=self.compression_opts)

        for  model_type, num_neurons_multiplier in models.items():
            start = time.time()
            model_str_def = build_model_string(model_type, num_neurons_multiplier)

            full_model_str_def = 'i%d_' % (self.ndays + self.nwindows) + model_str_def

            model_path_prefix = "mtl_%s" % (full_model_str_def)
            model, xscaler, yscaler = self.build_or_load_model(model_path_prefix, full_model_str_def, train_Y.shape)

            epochs = 50

            print("Model summary:")
            print(model.summary())

            if(xscaler is None or yscaler is None):
                print("Creating new scalers")

            xscaler, yscaler = annutils.create_or_update_xyscaler(xscaler, yscaler, train_X, train_Y)
            print("Xscaler Min[0]: %s" % xscaler.min_val[0])
            print("Xscaler Max[0]: %s" % xscaler.max_val[0])

            scaled_X = xscaler.transform(train_X)
            scaled_Y = yscaler.transform(train_Y)


            scaled_test_X = xscaler.transform(test_X)
            scaled_test_Y = yscaler.transform(test_Y)



            history = model.fit(
                scaled_X,
                scaled_Y,
                epochs=epochs,
                batch_size=36,
                validation_data=(scaled_test_X, scaled_test_Y),
                callbacks=[
                    keras.callbacks.EarlyStopping(
                        monitor="val_loss", patience=50, mode="min", restore_best_weights=True),
                    self.tensorboard_cb
                ],
                verbose=2,
            )

            # plot_history(history)

            model_savepath = os.path.join(local_root_path, "Experiments", self.experiment, 'models', model_path_prefix)
            # tf.saved_model.save(model, model_savepath)
            annutils.save_model(model_savepath, model, xscaler, yscaler)
            print('Model saved to %s' % model_savepath)
            print('Training time: %d min' % ((time.time() - start) / 60))

        print(f"Finished training model for {self.experiment}")

    def make_predictions(self, excel_files):
        print("Experiment: %s" % experiment)
        experiment_dir = os.path.join(local_root_path, "Experiments", experiment)

        model_dir = os.path.join(experiment_dir, "models")
        model_files = [f for f in os.listdir(model_dir) if f.endswith(".h5")]

        print("Local root: ",local_root_path)
        for data_file in tqdm(excel_files):
            print("Data file: %s" % data_file)
            data_path = os.path.join(local_root_path,data_file)
            dfinps, dfouts = annutils.read_and_split(data_path, self.num_sheets, self.observed_stations_ordered_by_median)
            for cn in dfinps.columns:
                print("Col "+cn)

            dfinps = annutils.create_antecedent_inputs(dfinps,ndays=self.ndays,window_size=self.window_size,nwindows=self.nwindows)
            dfinps, dfouts = annutils.synchronize(dfinps, dfouts)


            #get the name of the file without the extension
            # file_name = os.path.splitext(data_file)[0]
            file_name = os.path.splitext(os.path.basename(data_file))[0]

            dirs = ["input", "target", "prediction"]
            for dir in dirs:
                os.makedirs(os.path.join("Experiments", experiment, "results", dir), exist_ok=True)

            input_file = os.path.join("Experiments", experiment, "results", "input", file_name + ".csv")
            dfinps.to_csv(input_file, compression=self.compression_opts)

            # read_in = pd.read_csv(input_file, compression=compression_opts, index_col=0)

            target_file = os.path.join("Experiments", experiment, "results", "target", file_name + "_target.csv")
            dfouts.to_csv(target_file, compression=self.compression_opts)

            for model_file in tqdm(model_files):
                print("Model file: %s" % model_file)
                model_name = os.path.splitext(model_file)[0]

                model_prediction_dir = os.path.join(experiment_dir, "results", "prediction", model_name)
                os.makedirs(model_prediction_dir, exist_ok=True)

                location = os.path.join(model_dir, model_name)
                print("Location: %s" % location)
                #prediction = predict(location, dfinps, dfouts.columns)
                model=keras.models.load_model('%s.h5'%location,custom_objects={"mse_loss_masked": mse_loss_masked})
                xscaler,yscaler=joblib.load('%s_xyscaler.dump'%location)
                print("Xscaler Min[0]: %s" % xscaler.min_val[0])
                print("Xscaler Max[0]: %s" % xscaler.max_val[0])
                scaled_input = xscaler.transform(dfinps)
                print("scaled input")
                print(scaled_input.iloc[0])
                #print(scaled_input.iloc[1])
                #print(scaled_input)
                dfx = pd.DataFrame(scaled_input, dfinps.index, columns=dfinps.columns)

                yyp=model.predict(dfx, verbose=True)
                predicted_y = yscaler.inverse_transform(yyp)
                prediction = pd.DataFrame(predicted_y, index=dfinps.index, columns=dfouts.columns)
                print("prediction")
                print(prediction)
                prediction_file = os.path.join(model_prediction_dir, file_name + ".csv")
                print(f"Writing to {prediction_file}")
                prediction.to_csv(prediction_file, compression=self.compression_opts)

        print(f"Finished making predictions for {experiment}")

if __name__ == '__main__':

    from schimpy import schism_yaml
    import numpy as np
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    in_fname = "../../data/lathypcub_v1p1_ann_config.yaml"
    
    with open(in_fname, 'r') as f:
        # loader = RawLoader(stream)
        inputs = schism_yaml.load(f)
    dsp_home = inputs.get('dsp_home')
    if 'stas_include' in inputs.keys():
        stas_include = [v for v in inputs.get('stas_include')]
    else:
        stas_include = None
    
    # Initialize ###############################
    ann_mod = ModelANN(dsp_home, stas_include=stas_include)

    # Create Inputs ############################
    experiment = inputs.get('experiment') # name of folder where outputs are etc.
    
    if 'train_config' in inputs.keys():
        train_files = [tc.get('file').format(**globals(), **locals()) for tc in inputs.get('train_config')]
        train_windows = [[(tc.get('start'),tc.get('end'))] for tc in inputs.get('train_config')] # windows to use to extract training data from
    else:
        train_files = [t.format(**globals(), **locals()) for t in inputs.get('train_files')]
        train_windows = [[(s,e)] for s,e in zip(inputs.get('train_windows')[0].get('start'), 
                                              inputs.get('train_windows')[1].get('end'))] # windows to use to extract train data from
  
    if 'test_config' in inputs.keys():
        test_files = [tc.get('file').format(**globals(), **locals()) for tc in inputs.get('test_config')]
        test_windows = [[(tc.get('start'),tc.get('end'))] for tc in inputs.get('test_config')] # windows to use to extract training data from
    else:
        test_files = [t.format(**globals(), **locals()) for t in inputs.get('test_files')] # files to be used to test/train the model
        test_windows = [[(s,e)] for s,e in zip(inputs.get('test_windows')[0].get('start'), 
                                             inputs.get('test_windows')[1].get('end'))] # windows to use to extract test data from
        # fill any missing files or windows
        if len(test_files)==1 and len(test_windows)>1:
            test_files = np.repeat(test_files, len(test_windows))
        elif len(test_files)>1 and len(test_windows)==1:
            test_windows = np.repeat(test_windows, len(test_files))

    if 'vars_include' in inputs.keys():
        vars_include = [v for v in inputs.get('vars_include')]
    else:
        vars_include = None

    # run compile
    ann_mod.compile_inputs(experiment, train_files, train_windows, test_files, test_windows, 
                           vars_include=vars_include)

    # Train Model ##############################
    models = {m.get('name'):[p for p in m.get('params')] for m in inputs.get('models')}

    # run train
    ann_mod.train_models(models)

    # Make Predictions #########################
    excel_files = [e.format(**globals(),**locals()) for e in inputs.get('excel_files')] 
    
    # run predcit
    ann_mod.make_predictions(excel_files)
