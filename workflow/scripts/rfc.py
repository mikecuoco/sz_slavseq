#!/usr/bin/python3
__author__ = 'Ricardo S Jacomini'

import sys,gc
import traceback
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import itertools
import numpy as np
import pandas as pd
import re 

from pathlib import Path
import matplotlib.pyplot as plt

from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_score
from sklearn.metrics import f1_score
from sklearn.metrics import recall_score

from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder

import matplotlib.font_manager as fm
from matplotlib.collections import QuadMesh
import seaborn as sn

import warnings

warnings.filterwarnings("ignore")

sys.setrecursionlimit(10_000)

class TrainTest:
    def __init__(self, sample=None, exp="", input_type="pickle"):
        self.sample = snakemake.wildcards.donor # maybe change self.sample to self.donor
        self.dna_type = snakemake.wildcards.dna_type
        self.exp = exp
        self.input_type = input_type
        self.accurary_train = 0
        self.accurary_test = 0

    def train_classifier(self, features, labels):
        cla = RandomForestClassifier(bootstrap=True,
                                     n_estimators=100,
                                     oob_score=True,
                                     n_jobs=-1)

        cla.fit(features, labels)
        return cla

    def create_directory(self, directory):
        import os, errno

        try:
            os.makedirs(directory)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    def openDictionary_pickle(self, dict_data, infiles, printLoad=False):
        import csv
        import pandas as pd

        csv_columns = dict_data.keys()

        if printLoad == True: print("Loading :")

        for data in csv_columns:

            if printLoad == True: print("        " + data)
            
            # get appropriate split from input files
            files = [f for f in infiles if re.search(data, str(f))]

            # read split into dict_data
            dict_data[data] = pd.read_pickle(*files, compression='gzip')

        return dict_data

    def get_labelencode(self, df):

        le = LabelEncoder()

        classes = list(set(df))

        le.fit(classes)

        return le, le.fit_transform(df)

    def get_labelencode_inverse(self, le, df):

        return list(le.inverse_transform(df))

    def get_ypred(self, y_pred, proba, all_reads, le):

        df = y_pred.copy()

        df['all_reads_count'] = all_reads

        df.reset_index(inplace=True)

        for (a, b) in zip(le.classes_, proba):
            df[a + "_proba"] = proba[b]

        return df

    def get_groups(self, df):
        return df.groupby(['Y'], axis=0).size()

    def run_train_test(self, features, pred, phase, model=None, le=None, labels=None):

        if phase == "Training": model = self.train_classifier(features, labels)

        predictions = model.predict(features)

        df = accuracy_score(labels, predictions, sample_weight=None)

        proba = model.predict_proba(features)

        pred['Y_pred'] = self.get_labelencode_inverse(le, predictions)

        pred = self.get_ypred(pred, pd.DataFrame(proba), features['all_reads.count'], le)

        log_stdout = open(self.outDir + "/%s_report_performance.txt" % phase, "w+")
        print("\n Report performance %s" % phase)
        self.report_performance(model, features, labels, predictions, log_stdout, phase, le)

        log_stdout = open(self.outDir + "/%s_report_feature_importance.txt" % phase, "w+")
        print("\n Report feature importance")
        self.report_feature_importances(model, features.columns, log_stdout)

        return df, pred, model

    def run(self):

        print("Sample %s" % self.sample) # maybe donor wildcard?

        for k in range(snakemake.params.num_folds):

            # TODO make this robust to output being a directory or a file
            for f in snakemake.output:
                if re.search(f'fold_{k}', f):
                    self.outDir = str(Path(f).parent)
            
            self.inFiles = [Path(f) for f in snakemake.input if re.search(f'fold_{k}', f)]
            
            # TODO: examine why this if statement is not entered; why does the path exist?
            # if not Path(self.outDir).exists():
            self.create_directory(self.outDir)

            # clear_output()

            # Loading dataset ------------------------------------------------------------------------------------

            train_x, train_y, test_x, test_y = self.get_sample(self.inFiles)

            le, labels = self.get_labelencode(train_y['Y'])

            # classes
            self.classes = list(le.classes_)

            if self.get_groups(train_y).count() != self.get_groups(test_y).count():
                log_stdout = open(self.outDir + "/" + self.sample + "_has_no_equal_classes.txt", "w+")
                print("\nThis cell (%s) fold (%i) has no enough classes to classifier... " % (self.sample, k),
                        file=log_stdout, flush=True)
                print("\nTrain classes (%d) - Test classes (%d)" % (
                self.get_groups(train_y).count(), self.get_groups(test_y).count()), file=log_stdout, flush=True)
                continue

            # Training the model -----------------------------------------------------------------------------

            print("\nTraining the model... ")
            Train, df, trained_model = self.run_train_test(train_x.copy(), train_y.copy(),
                                                            phase='Training',
                                                            le=le,
                                                            labels=labels)

            df.to_csv(self.outDir + "/Training_y_pred.csv", index=False, header=True)

            # Testing the model ----------------------------------------------------------------------------------

            print("\nTesting the model... ")

            le, labels = self.get_labelencode(test_y['Y'])
            Test, df, trained_model = self.run_train_test(test_x.copy(), test_y.copy(), phase='Testing',
                                                            model=trained_model, le=le, labels=labels)

            df.to_csv(self.outDir + "/Testing_y_pred.csv", index=False, header=True)

            # report sys.out --- 
            print("\nTrain Accuracy :: {} - Test Accuracy  :: {}".format(Train, Test))

            self.accurary_train = Train
            self.accurary_test = Test

            # saving ------------------------------------------------------------------------------------------ 
            print("\nSaving...")

            df = pd.DataFrame(data=None, columns={'Train', 'Test'}, dtype='float64')

            df = df.append([{'Train': Train, 'Test': Test}])
            df.to_csv(self.outDir + "/Train_Test_Accuracy.csv", sep=';', index=False, header=True)

        print("%s : finishing" % (self.sample))

    def get_sample(self, files):

        print("Processing cell ({})... \n".format(self.sample))

        print("Loading dataset...")

        df = {"X_train": None, "X_test": None, "Y_train": None, "Y_test": None}

        TrainTest = self.openDictionary_pickle(df, files)
        train_y = TrainTest["Y_train"].to_frame('Y')
        test_y = TrainTest["Y_test"].to_frame('Y')

        train_x = TrainTest["X_train"]
        test_x = TrainTest["X_test"]

        # Train and Test dataset size details
        print("Train_x Shape :: ", train_x.shape)
        print("Train_y Shape :: ", train_y.shape)
        print("Test_x Shape :: ", test_x.shape)
        print("Test_y Shape :: ", test_y.shape)
        print("\n")

        return train_x, train_y, test_x, test_y

    def __str__(self):
        return "Sample %s - Tain Accuracy: %s - Test Accuracy %s " % (
        self.sample, self.accurary_train, self.accurary_test)

    def __repr__(self):
        return str(self.__dict__)

    def __hash__(self):
        return hash(self.sample)

    @classmethod
    def get_classname(cls):
        return cls.__name__

    def use_classname(self):
        return self.get_classname()

    def report_feature_importances(self, cla, feature_names, log_stdout):
        ii = np.argsort(cla.feature_importances_)
        print("Feature importance:", file=log_stdout, flush=True)
        for i in range(len(ii)):
            print('', '\t'.join([str(cla.feature_importances_[ii[i]]), feature_names[ii[i]]]), file=log_stdout,
                  flush=True)
        print(file=log_stdout, flush=True)

    def report_performance(self, cla, features, labels, pred, log_stdout, title, le):

        class_names = cla.classes_
        cnf_matrix = confusion_matrix(labels, pred)

        acc = np.trace(cnf_matrix) / float(np.sum(cnf_matrix))
        pred_val_axis = 'x'
        cmap = 'Oranges'
        fz = 10
        figsize = [5, 5]

        show_null_values = 2
        # cols_names = ["RL1","KNRGL","OTHER"]
        cols_names = list(le.classes_)

        df_cm = pd.DataFrame(cnf_matrix, index=cols_names, columns=cols_names)

        self.pretty_plot_confusion_matrix(df_cm, accuracy=acc, fz=fz, cmap=cmap, figsize=figsize,
                                          show_null_values=show_null_values,
                                          pred_val_axis=pred_val_axis,
                                          title=('Confusion matrix (%s) ' % title),
                                          fn=self.outDir + "/" + title + '_Confusion-pretty.png')

        self.plot_confusion_matrix(cnf_matrix,
                                   title=('Confusion matrix (%s), without normalization' % title),
                                   fn=self.outDir + "/" + title + '_Confusion.png')

        self.plot_confusion_matrix(cnf_matrix, normalize=True,
                                   title=('Normalized confusion matrix (%s)' % title),
                                   fn=self.outDir + "/" + title + '_Confusion-Normalized.png')

        print("Classifier:", file=log_stdout, flush=True)
        print(cla, file=log_stdout, flush=True)
        print(file=log_stdout, flush=True)
        print(("%s set confusion matrix:" % title), file=log_stdout, flush=True)

        print(cnf_matrix, file=log_stdout, flush=True)
        print(file=log_stdout, flush=True)

        # https://en.wikipedia.org/wiki/Precision_and_recall
        # https://machinelearningmastery.com/roc-curves-and-precision-recall-curves-for-classification-in-python/
        print("Cross-validation metrics:", file=log_stdout, flush=True)
        for jj in ['precision', 'recall', 'f1']:

            if jj == 'precision':
                scores = precision_score(labels, pred, average='weighted')
            elif jj == 'recall':
                scores = recall_score(labels, pred, average='weighted')
            else:
                scores = f1_score(labels, pred, average='weighted')

            print('', jj, scores, file=log_stdout, flush=True)

            print('\t', 'mean', jj, scores.mean(), file=log_stdout, flush=True)
            print('\t', 'std', jj, scores.std(), file=log_stdout, flush=True)

        print(file=log_stdout, flush=True)

    def plot_confusion_matrix(self, cm,
                              normalize=False,
                              title='Confusion matrix',
                              cmap=plt.cm.Blues, fn="cm"):
        """
        This function prints and plots the confusion matrix.
        Normalization can be applied by setting `normalize=True`.
        """
        if normalize:
            cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
            print("Normalized confusion matrix")
        else:
            print('Confusion matrix, without normalization')

        print(cm)

        plt.imshow(cm, interpolation='nearest', cmap=cmap)
        plt.title(title)
        plt.colorbar()
        tick_marks = np.arange(len(self.classes))
        plt.xticks(tick_marks, self.classes, rotation=45)
        plt.yticks(tick_marks, self.classes)

        fmt = '.2f' if normalize else 'd'
        thresh = cm.max() / 2.
        for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
            plt.text(j, i, format(cm[i, j], fmt),
                     horizontalalignment="center",
                     color="white" if cm[i, j] > thresh else "black")

        plt.tight_layout()
        plt.ylabel('True label')
        plt.xlabel('Predicted label')

        plt.savefig(fn, bbox_inches='tight')
        plt.close()

    def configcell_text_and_colors(self, array_df, lin, col, oText, facecolors, posi, fz, fmt, show_null_values=0):
        """
          config cell text and colors
          and return text elements to add and to dell
          @TODO: use fmt
        """
        text_add = []
        text_del = []
        cell_val = array_df[lin][col]
        tot_all = array_df[-1][-1]
        per = (float(cell_val) / tot_all) * 100
        curr_column = array_df[:, col]
        ccl = len(curr_column)

        # last line  and/or last column
        if (col == (ccl - 1)) or (lin == (ccl - 1)):
            # tots and percents
            if cell_val != 0:
                if (col == ccl - 1) and (lin == ccl - 1):
                    tot_rig = 0
                    for i in range(array_df.shape[0] - 1):
                        tot_rig += array_df[i][i]
                    per_ok = (float(tot_rig) / cell_val) * 100
                elif col == ccl - 1:
                    tot_rig = array_df[lin][lin]
                    per_ok = (float(tot_rig) / cell_val) * 100
                elif lin == ccl - 1:
                    tot_rig = array_df[col][col]
                    per_ok = (float(tot_rig) / cell_val) * 100
                per_err = 100 - per_ok
            else:
                per_ok = per_err = 0

            per_ok_s = ['%.2f%%' % per_ok, '100%'][per_ok == 100]

            # text to DEL
            text_del.append(oText)

            # text to ADD
            font_prop = fm.FontProperties(weight='bold', size=fz)
            text_kwargs = dict(color='w', ha="center", va="center", gid='sum', fontproperties=font_prop)
            lis_txt = ['%d' % cell_val, per_ok_s, '%.2f%%' % per_err]
            lis_kwa = [text_kwargs]
            dic = text_kwargs.copy()
            dic['color'] = 'g'
            lis_kwa.append(dic)
            dic = text_kwargs.copy()
            dic['color'] = 'r'
            lis_kwa.append(dic)
            lis_pos = [(oText._x, oText._y - 0.3), (oText._x, oText._y), (oText._x, oText._y + 0.3)]
            for i in range(len(lis_txt)):
                newText = dict(x=lis_pos[i][0], y=lis_pos[i][1], text=lis_txt[i], kw=lis_kwa[i])
                # print 'lin: %s, col: %s, newText: %s' %(lin, col, newText)
                text_add.append(newText)
            # print '\n'

            # set background color for sum cells (last line and last column)
            carr = [0.27, 0.30, 0.27, 1.0]
            if (col == ccl - 1) and (lin == ccl - 1):
                carr = [0.17, 0.20, 0.17, 1.0]
            facecolors[posi] = carr

        else:
            if per > 0:
                txt = '%s\n%.2f%%' % (cell_val, per)
            else:
                if show_null_values == 0:
                    txt = ''
                elif show_null_values == 1:
                    txt = '0'
                else:
                    txt = '0\n0.0%'
            oText.set_text(txt)

            # main diagonal
            if col == lin:
                # set color of the text in the diagonal to white
                oText.set_color('w')
                # set background color in the diagonal to blue
                facecolors[posi] = [0.35, 0.8, 0.55, 1.0]
            else:
                oText.set_color('r')

        return text_add, text_del


    def insert_totals(self, df_cm):
        """ insert total column and line (the last ones) """
        sum_col = []
        for c in df_cm.columns:
            sum_col.append(df_cm[c].sum())
        sum_lin = []
        for item_line in df_cm.iterrows():
            sum_lin.append(item_line[1].sum())
        df_cm['total_rows'] = sum_lin
        sum_col.append(np.sum(sum_lin))
        df_cm.loc['total_cols'] = sum_col
        # print ('\ndf_cm:\n', df_cm, '\n\b\n')

    def get_new_fig(self, fn, figsize=[9, 9]):
        """ Init graphics """
        fig1 = plt.figure(fn, figsize)
        ax1 = fig1.gca()  # Get Current Axis
        ax1.cla()  # clear existing plot
        return fig1, ax1

    #
    # https://github.com/wcipriano/pretty-print-confusion-matrix/blob/master/confusion_matrix_pretty_print.py
    def pretty_plot_confusion_matrix(self, df_cm, accuracy=None, annot=True, cmap="Oranges", fmt='.2f', fz=11,
                                     lw=0.5, cbar=False, figsize=[8, 8], show_null_values=0, pred_val_axis='y',
                                     title='Confusion matrix', fn="cm"):
        """
          print conf matrix with default layout (like matlab)
          params:
          df_cm          dataframe (pandas) without totals
          annot          print text in each cell
          cmap           Oranges,Oranges_r,YlGnBu,Blues,RdBu, ... see:
          fz             fontsize
          lw             linewidth
          pred_val_axis  where to show the prediction values (x or y axis)
                         'col' or 'x': show predicted values in columns (x axis) instead lines
                         'lin' or 'y': show predicted values in lines   (y axis)
        """
        if pred_val_axis in ('col', 'x'):
            xlbl = 'Predicted label'
            ylbl = 'Actual label'
        else:
            xlbl = 'Actual label'
            ylbl = 'Predicted label'
            df_cm = df_cm.T

        # create "Total" column
        self.insert_totals(df_cm)

        if cmap is None:
            cmap = plt.get_cmap('Blues')

        # this is for print always in the same window
        fig, ax1 = self.get_new_fig('Conf matrix default', figsize)

        # thanks for seaborn
        ax = sn.heatmap(df_cm, annot=annot, annot_kws={"size": fz}, linewidths=lw, ax=ax1,
                        cbar=cbar, cmap=cmap, linecolor='w', fmt=fmt)

        # set ticklabels rotation
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, fontsize=10)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=25, fontsize=10)

        # Turn off all the ticks
        for t in ax.xaxis.get_major_ticks():
            t.tick1On = False
            t.tick2On = False
        for t in ax.yaxis.get_major_ticks():
            t.tick1On = False
            t.tick2On = False

        # face colors list
        quadmesh = ax.findobj(QuadMesh)[0]
        facecolors = quadmesh.get_facecolors()

        # iter in text elements
        array_df = np.array(df_cm.to_records(index=False).tolist())
        text_add = []
        text_del = []
        posi = -1  # from left to right, bottom to top.
        for t in ax.collections[0].axes.texts:  # ax.texts:
            pos = np.array(t.get_position()) - [0.5, 0.5]
            lin = int(pos[1])
            col = int(pos[0])
            posi += 1
            # print ('>>> pos: %s, posi: %s, val: %s, txt: %s' %(pos, posi, array_df[lin][col], t.get_text()))

            # set text
            txt_res = self.configcell_text_and_colors(array_df, lin, col, t, facecolors, posi, fz, fmt,
                                                      show_null_values)

            text_add.extend(txt_res[0])
            text_del.extend(txt_res[1])

        # remove the old ones
        for item in text_del:
            item.remove()
        # append the new ones
        for item in text_add:
            ax.text(item['x'], item['y'], item['text'], **item['kw'])

        # titles and legends
        ax.set_title(title)
        ax.set_xlabel(xlbl)
        ax.set_ylabel(ylbl)

        plt.xlabel('Predicted label\nAccuracy={:0.4f}; Misclass={:0.4f}'.format(accuracy, 1 - accuracy))
        plt.tight_layout()  # set layout slim

        plt.savefig(fn, bbox_inches='tight')
        plt.close()       

    def plot_confusion_matrix(self, cm,
                              normalize=False,
                              title='Confusion matrix',
                              cmap=plt.cm.Blues, fn="cm"):
        """
        This function prints and plots the confusion matrix.
        Normalization can be applied by setting `normalize=True`.
        """
        if normalize:
            cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
            print("Normalized confusion matrix")
        else:
            print('Confusion matrix, without normalization')

        print(cm)

        plt.imshow(cm, interpolation='nearest', cmap=cmap)
        plt.title(title)
        plt.colorbar()
        tick_marks = np.arange(len(self.classes))
        plt.xticks(tick_marks, self.classes, rotation=45)
        plt.yticks(tick_marks, self.classes)

        fmt = '.2f' if normalize else 'd'
        thresh = cm.max() / 2.
        for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
            plt.text(j, i, format(cm[i, j], fmt),
                     horizontalalignment="center",
                     color="white" if cm[i, j] > thresh else "black")

        plt.tight_layout()
        plt.ylabel('True label')
        plt.xlabel('Predicted label')

        plt.savefig(fn, bbox_inches='tight')
        plt.close()

if __name__ == '__main__':
 
    try:
        train_test = TrainTest()
        train_test.run()

    except:  # catch *all* exceptions
        sys.stderr = open(snakemake.log[0], 'w')
        traceback.print_exc()
        sys.stderr.close()

    finally:
        # cleanup code in here
        gc.collect()

    sys.stdout.close()
