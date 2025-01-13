import PySimpleGUI as sg
import numpy as np
import pandas as pd
from tkinter import messagebox, simpledialog
from datetime import datetime
import random
import copy
import CCMPS_actions


## parameters:
#------------------------------------------------------------------------------




## internal functions:
#------------------------------------------------------------------------------

def fract_refreshtable (window, table):
    window['-fractionation_table-'].Update(values = table)
    return

def tp_refreshtable (window, table):
    window['-tp_table-'].Update(values = table)
    return


def fract_modifytable (title, prompt, values, fract_tables, pos, q, ask):
    if values['-fractionation_table-']:
        path = values['-fractionation_path-']
        table = fract_tables[path]
        if ask == 'integer':
            value = simpledialog.askinteger(title, prompt)
            p = 0
            if value:
                for i in values['-fractionation_table-']:
                    table[i][pos] = value +p
                    p = p+q
                fract_tables[path] = table
        elif ask == 'string':
            value = simpledialog.askstring(title, prompt)
            if value:
                for i in values['-fractionation_table-']:
                    table[i][pos] = value
                fract_tables[path] = table
    else:
        messagebox.showerror("Error", "Select (a) sample(s).")
    return(values, fract_tables)


def fract_buttons (window, status):
    active = ['-fractionation_add-',
              '-fractionation_remove-',
              '-fractionation_edit_remove-',
              '-fractionation_edit_keep-',
              '-fractionation_edit_condition-',
              '-fractionation_edit_replicate-',
              '-fractionation_edit_fractions-',
              '-fractionation_edit_identifier-',
              '-fractionation_parameters-',
              '-fractionation_start-']
    inactive = ['-fractionation_reset-',
                '-fractionation_summary-',
                ]
    for button in active:
        window[button].Update(disabled = status)
    for button in inactive:
        window[button].Update(disabled = not status)
    if status:
        window['-fractionation_status-'].Update(value = 'done!', text_color = 'dark green')
    else:
        window['-fractionation_status-'].Update(value = '...ready!', text_color = 'white')
    return


def tp_buttons (window, status):
    active = ['-tp_add-',
              '-tp_remove-',
              '-tp_edit_remove-',
              '-tp_edit_keep-',
              '-tp_edit_condition-',
              '-tp_edit_identifier-',
              '-tp_start-',
              '-tp_parameters-']
    inactive = ['-tp_reset-',
                '-tp_summary-',
                '-tp_export-']
    for button in active:
        window[button].Update(disabled = status)
    for button in inactive:
        window[button].Update(disabled = not status)
    if status:
        window['-tp_status-'].Update(value = 'done!', text_color = 'dark green')
    else:
        window['-tp_status-'].Update(value = '...ready!', text_color = 'white')

def reset_fract ():
    data_ways = {'class' : [],
                 'vis' : []}
    std_ways = {'class' : [],
                'vis' : []}
    intermediate_data = {}
    protein_info = {}
    conditions = []
    return data_ways, std_ways, intermediate_data, protein_info, conditions


def reset_tp():
    tp_data = {}
    tp_intermediate = {}
    tp_info = pd.DataFrame()
    tp_conditions = []
    tp_icorr = {}
    return tp_data, tp_intermediate, tp_info, tp_conditions, tp_icorr


def resetinput ():
    input_paths = []
    input_tables = {}
    input_data = {}
    ident_pos = {}
    return input_paths, input_tables, input_data, ident_pos


def fract_clearinput (window):
    window['-fractionation_path-'].Update(values = [])
    window['-fractionation_table-'].Update(values = [])
    return

def tp_clearinput (window):
    window['-tp_path-'].Update (values = [])
    window['-tp_table-'].Update (values = [])
    return




## external functions:
#------------------------------------------------------------------------------

def reset_infract ():
    fract_indata = {}
    fract_identifiers = {}
    return fract_indata, fract_identifiers



def reset_intp ():
    tp_indata = {}
    tp_identifiers = {}
    return tp_indata, tp_identifiers


def fract_add (values, window, fract_paths, fract_tables, fract_indata, fract_pos, fract_identifiers):
    filename = sg.popup_get_file('Chose dataset', no_window = True, file_types=(('Tab Separated Values', '*.tsv'),('Text (tab delimited)', '*.txt')))
    if filename:
        fract_paths.append(filename)
        window['-fractionation_path-'].Update(values = fract_paths, value = filename)
        data = pd.read_csv(filename, sep = "\t", header = 0)
        data = data.replace('NaN', np.nan)
        data = data.replace('Filtered', np.nan)
        colnames = data.columns.values.tolist()
        table = []
        for name in colnames:
            namelist = [name,'','','']
            table.append(namelist)
        fract_tables[filename] = table
        fract_indata[filename] = data
        fract_pos[filename] = []
        fract_identifiers[filename] = []
        
        fract_refreshtable(window, table)
    return


def fract_rem (values, window, fract_paths, fract_tables, fract_data):
    sure = sg.popup_yes_no('Remove data from list?')
    if sure == 'Yes':
        fract_paths.remove(values['-fractionation_path-'])
        del fract_tables[values['-fractionation_path-']]
        if fract_paths:
            curr = fract_paths[0]
            fract_refreshtable(window, fract_tables[curr])
        else:
            curr = []
            fract_refreshtable(window, curr)
        window['-fractionation_path-'].Update(values = fract_paths, value = curr)
    else:
        pass
    return


def fract_defrem (values, window, fract_tables):
    path = values['-fractionation_path-']
    selected = values['-fractionation_table-']
    table = fract_tables[path]
    for index in sorted(selected, reverse = True):
        del table[index]
    fract_tables[path] = table
    window['-fractionation_table-'].Update(values = fract_tables[path])
    return


def fract_defkeep (values, window, fract_tables):
    path = values['-fractionation_path-']
    table = fract_tables[path]
    for pos in values['-fractionation_table-']:
        table[pos][1] = '[KEEP]'
        table[pos][2] = '-'
        table[pos][3] = '-'
    fract_tables[path] = table
    window['-fractionation_table-'].Update(values = fract_tables[path])
    return


def fract_defcon (values, window, fract_tables):
    values, fract_tables = fract_modifytable('Set Condition', 'Condition Name:', values, fract_tables, 1, 0, 'string')
    window['-fractionation_table-'].Update(values = fract_tables[values['-fractionation_path-']])
    return


def fract_defrep (values, window, fract_tables):
    values, fract_tables = fract_modifytable('Set Replicate', 'Replicate Number:', values, fract_tables, 2, 0, 'integer')
    window['-fractionation_table-'].Update(values = fract_tables[values['-fractionation_path-']])
    return


def fract_deffract (values, window, fract_tables):
    values, fract_tables = fract_modifytable('Set Fractions', 'FIRST Fraction Number:', values, fract_tables, 3, 1, 'integer')
    window['-fractionation_table-'].Update(values = fract_tables[values['-fractionation_path-']])
    return


def fract_defident (values, window, input_tables, ident_pos, identifiers):
    pos = values['-fractionation_table-']
    if pos:
        if len(pos) > 1:
            messagebox.showerror("Error", "Please set only one Identifier!")
        elif len(pos) == 1:
            path = values['-fractionation_path-']
            table = input_tables[path]
            if ident_pos[path]:
                table[ident_pos[path][0]][1] = ''
                table[ident_pos[path][0]][2] = ''
                table[ident_pos[path][0]][3] = ''
            identifiers[path] = table[pos[0]][0]
            ident_pos[path] = pos
            table[pos[0]][1] = '[IDENTIFIER]'
            table[pos[0]][2] = '-'
            table[pos[0]][3] = '-'
            input_tables[path] = table
            window['-fractionation_table-'].Update(values = input_tables[values['-fractionation_path-']])
        else:
            messagebox.showerror("Error", "No sample selected.")
    return identifiers


def fract_export (values, data, protein_info):
    export_folder = sg.popup_get_folder('Export Folder')
    if export_folder:
        experiment = simpledialog.askstring('Export', 'Experiment Name: ')
        now = datetime.now()
        time = now.strftime("%Y%m%d%H%M%S")
        # export_full = {'class' : pd.DataFrame(),
        #                'vis' : pd.DataFrame()}
        for way in data:
            data_way = pd.DataFrame()
            for condition in data[way]:
                path = export_folder + '/' + time + '_' + experiment + '_' + way + '_' + condition + '.txt'
                export_data = data[way][condition]
                for info in protein_info:
                    export_data = pd.merge(export_data, protein_info[info], left_index = True, right_index = True, how = 'left')
                export_data.to_csv(path, header = True, index = True, index_label = 'Identifier', sep = '\t', mode = 'a')
                
                data_way = pd.merge(data_way, data[way][condition], left_index = True, right_index = True, how = 'outer')
                
                
                # for replicate in data[way][condition]:
                #     path = export_folder + '/' + time + '_' + experiment + '_' + way + '_' + condition + '_' + replicate + '.txt'
                #     data[way][condition][replicate].to_csv(path, header = True, index = True, index_label = 'Identifier', sep = '\t', mode = 'a')
                #     export_full[way] = pd.merge(export_full[way], data[way][condition][replicate], left_index = True, right_index = True, how = 'outer')
            
            # for info in protein_info:
            #     export_full[way] = pd.merge(export_full[way], protein_info[info], left_index = True, right_index = True, how = 'left')
            
            for info in protein_info:
                data_way = pd.merge(data_way, protein_info[info], left_index = True, right_index = True, how = 'left')
            
            
            path = export_folder + '/' + time + '_' + experiment + '_' + way + '_' + 'COMBINED' + '.txt'
            data_way.to_csv(path, header = True, index = True, index_label = 'Identifier', sep = '\t', mode = 'a')
    return




#------------------------------------------------------------------------------

def is_float (element):
    try:
        float(element)
        return True
    except:
        return False

def convert_to_float(x):
    try:
        return float(x)
    except ValueError:
        return x

def tp_add (values, window, tp_paths, tp_tables, tp_indata, tp_pos, tp_identifiers):
    filename = sg.popup_get_file('Chose dataset', no_window = True, file_types=(('Tab Separated Values', '*.tsv'),('Text (tab delimited)', '*.txt')))
    if filename:
        tp_paths.append(filename)
        window['-tp_path-'].Update(values = tp_paths, value = filename)
        data = pd.read_csv(filename, sep = "\t", header = 0)
        data = data.replace('NaN', np.nan)
        data = data.replace('Filtered', np.nan)
        data = data.applymap(convert_to_float)
        
        rows_with_float = data.applymap(is_float).any(axis = 1)
        data = data[rows_with_float]
        
        colnames = data.columns.values.tolist()
        table = []
        for name in colnames:
            namelist = [name, '']
            table.append(namelist)
        tp_tables[filename] = table
        tp_indata[filename] = data
        tp_pos[filename] = []
        tp_identifiers[filename] = []
        tp_refreshtable(window, table)
    return




def tp_rem (values, window, tp_paths, tp_tables, tp_data):
    sure = sg.popup_yes_no('Remove data from list?')
    if sure == 'Yes':
        tp_paths.remove(values['-tp_path-'])
        del tp_tables[values['-tp_path-']]
        #del tp_data[values['-tp_path-']]
        if tp_paths:
            curr = tp_paths[0]
            tp_refreshtable(window, tp_tables[curr])
        else:
            curr = []
            tp_refreshtable(window, curr)
        window['-tp_path-'].Update(values = tp_paths, value = curr)
    else:
        pass
    return


def tp_defrem (values, window, tp_tables):
    path = values['-tp_path-']
    selected = values['-tp_table-']
    table = tp_tables[path]
    for index in sorted(selected, reverse = True):
        del table[index]
    tp_tables[path] = table
    window['-tp_table-'].Update(values = tp_tables[path])
    return


def tp_defkeep (values, window, tp_tables):
    path = values['-tp_path-']
    table = tp_tables[path]
    for pos in values['-tp_table-']:
        table[pos][1] = '[KEEP]'
    tp_tables[path] = table
    window['-tp_table-'].Update(values = tp_tables[path])
    return


def tp_defcon (values, window, tp_tables):
    if values['-tp_table-']:
        path = values['-tp_path-']
        table = tp_tables[path]
        value = simpledialog.askstring('Set Condition', 'Condition Name')
        if value:
            for i in values['-tp_table-']:
                table[i][1] = value
            tp_tables[path] = table
    else:
        messagebox.showerror("Error", "Select (a) sample(s).")
    window['-tp_table-'].Update(values = tp_tables[values['-tp_path-']])
    return


def tp_defident (values, window, tp_tables, tp_pos, tp_identifiers):
    pos = values['-tp_table-']
    if pos:
        if len(pos) > 1:
            messagebox.showerror("Error", "Please set only one Identifier!")
        elif len(pos) == 1:
            path = values['-tp_path-']
            table = tp_tables[path]
            if tp_pos[path]:
                table[tp_pos[path][0]][1] = ''
            tp_identifiers[path] = table[pos[0]][0]
            tp_pos[path] = pos
            table[pos[0]][1] = '[IDENTIFIER]'
            tp_tables[path] = table
            window['-tp_table-'].Update(values = tp_tables[values['-tp_path-']])
        else:
            messagebox.showerror("Error", "No sample selected.")
    return tp_identifiers


def tp_export (tp_data, tp_info):
    export_folder = sg.popup_get_folder('Export Folder')
    if export_folder:
        experiment = simpledialog.askstring('Export', 'Experiment Name: ')
        now = datetime.now()
        time = now.strftime("%Y%m%d%H%M%S")
        export_full = pd.DataFrame()
        for condition in tp_data:
            path = export_folder + '/' + time + '_' + experiment + '_' + condition + '.txt'
            tp_data[condition].to_csv(path, header = True, index = True, index_label = 'Identifier', sep = '\t', mode = 'a')
            export_full = pd.merge(export_full, tp_data[condition], left_index = True, right_index = True, how = 'outer')
        for info in tp_info:
            export_full = pd.merge(export_full, tp_info[info], left_index = True, right_index = True, how = 'left')
        path = export_folder + '/' + time + '_' + experiment + '_' + 'combined' + '.txt'
        export_full.to_csv(path, header = True, index = True, index_label = 'Identifier', sep = '\t', mode = 'a')
    return








#------------------------------------------------------------------------------




def check_markers(marker_sets):
    is_markers = True
    
    if marker_sets:
        for file in marker_sets:
            if marker_sets[file]['identifier_col'] == '-' or marker_sets[file]['class_col'] == '-':
                is_markers = False
    else:
        is_markers = False
    
    
    
    
    # print(marker_sets[file])
    return(is_markers)




def refresh_markertable (window, values, marker_sets):
    file_list = []
    for markerfile in marker_sets:
        file_list.append(markerfile)
    window['-marker_list-'].Update(values = file_list)
    if file_list:
        window['-marker_list-'].Update(set_to_index = 0)
        window['-marker_key-'].Update(values = marker_sets[file_list[0]]['table'].columns.tolist(), value = marker_sets[file_list[0]]['identifier_col'])
        window['-marker_class-'].Update(values = marker_sets[file_list[0]]['table'].columns.tolist(), value = marker_sets[file_list[0]]['class_col'])
    return


def refresh_markercols (window, values, marker_sets):
    try:
        window['-marker_key-'].Update(values = marker_sets[values['-marker_list-'][0]]['table'].columns.tolist(),
                                                 value = marker_sets[values['-marker_list-'][0]]['identifier_col'])
        window['-marker_class-'].Update(values = marker_sets[values['-marker_list-'][0]]['table'].columns.tolist(),
                                                 value = marker_sets[values['-marker_list-'][0]]['class_col'])
    except:
        window['-marker_key-'].Update(values = [],
                                      value = '-')
        window['-marker_class-'].Update(values = [],
                                      value = '-')
    return









def marker_add (window, values, marker_sets):
    filename = sg.popup_get_file('Select a new Marker List!', no_window = True, file_types=(('Tab delimited Text', '*.txt'),('Tab Separated Values', '*.tsv')))
    if filename:
        marker_sets[filename] = {}
        # marker_sets[filename]['table'] = pd.read_csv(filename, sep = "\t", header = 0).apply(lambda x: x.astype(str))
        marker_sets[filename]['table'] = pd.read_csv(filename, sep = "\t", header = 0).apply(lambda x: x.astype(str).str.upper())
        marker_sets[filename]['identifier_col'] = '-'
        marker_sets[filename]['class_col'] = '-'
        marker_sets[filename]['classes'] = []
        refresh_markertable(window, values, marker_sets)
        
        # window['-marker_test-'].Update(disabled = False)
        # window['-marker_profiles-'].Update(disabled = False)
        # window['-marker_remove-'].Update(disabled = False)
    return






def marker_remove (window, values, marker_sets):
    del marker_sets[values['-marker_list-'][0]]
    refresh_markertable(window, values, marker_sets)
    if not len(marker_sets) > 0:
        window['-marker_test-'].Update(disabled = True)
        window['-marker_profiles-'].Update(disabled = True)
        window['-marker_remove-'].Update(disabled = True)
    return marker_sets




def marker_setkey (values, marker_sets):
    marker_sets[values['-marker_list-'][0]]['identifier_col'] = values['-marker_key-']
    return


def marker_setclass (values, marker_sets):
    marker_sets[values['-marker_list-'][0]]['class_col'] = values['-marker_class-']
    marker_sets[values['-marker_list-'][0]]['classes'] = list(set(marker_sets[values['-marker_list-'][0]]['table'][values['-marker_class-']]))
    marker_conv = create_conversion(marker_sets)
    return marker_conv



def create_conversion(marker_sets):
    marker_conv = {}
    for path in marker_sets:
        for classname in marker_sets[path]['classes']:
            marker_conv[classname] = classname
    return marker_conv



def create_markerlist (marker_sets, marker_conv, marker_params):
    markerset = pd.DataFrame(columns = ['name'])
    counter = 1
    for path in marker_sets:
        mset = marker_sets[path]['table'][[marker_sets[path]['identifier_col'], marker_sets[path]['class_col']]]
        mset.rename(columns = {marker_sets[path]['identifier_col'] : 'name', marker_sets[path]['class_col'] : 'class' + str(counter)}, inplace = True)
        for classname in marker_conv:
            mset['class' + str(counter)].replace({classname : marker_conv[classname]}, inplace = True)
            mset['class' + str(counter)].replace(r'^\s*$', np.nan, regex=True, inplace = True)
            mset = mset[mset['class' + str(counter)].notna()]
        markerset = pd.merge(markerset, mset, on = 'name', how = 'outer')
        counter +=1
    markerset.set_index('name', inplace = True)
    if marker_params['what'] == 'unite':
        pass
    elif marker_params['what'] == 'intersect':
        markerset.dropna(inplace = True)
    if marker_params['how'] == 'majority':
        markerset_final = pd.DataFrame(markerset.mode(axis = 1, dropna = True)[0]).rename(columns = {0 : 'class'})
    elif marker_params['how'] == 'exclude':
        markerset_final = markerset.mode(axis = 1, dropna = True).fillna(np.nan)
        if 1 in markerset_final.columns:
            markerset_final = pd.DataFrame(markerset_final[markerset_final[1].isnull()][0]).rename(columns = {0 : 'class'})
        else:
            markerset_final.rename(columns = {0 : 'class'}, inplace = True)
    return markerset_final

#------------------------------------------------------------------------------

def session_save (fract_paths, fract_tables, fract_pos, fract_indata, fract_data, fract_std, fract_intermediate, fract_identifiers, fract_info, fract_preparams, tp_paths, tp_tables, tp_pos, tp_indata, tp_data, tp_intermediate, tp_identifiers, tp_info, tp_preparams, marker_sets, marker_params, marker_conv, fract_conditions, fract_full, fract_full_up, fract_marker, fract_marker_vis, fract_marker_up, fract_mixed_up, fract_test, svm_marker, svm_test, svm_metrics, marker_list, learning_xyz, NN_params, results, comparison, marker_fractkey, status):
    filename = sg.popup_get_file('Save Session', no_window = True, file_types=(('Numpy', '*.npy'),), save_as = True)
    if filename:
        file = {'fract_paths' : fract_paths,
                'fract_tables' : fract_tables,
                'fract_pos' : fract_pos,
                'fract_indata' : fract_indata,
                'fract_data' : fract_data,
                'fract_std' : fract_std,
                'fract_intermediate' : fract_intermediate,
                'fract_identifiers' : fract_identifiers,
                'fract_info' : fract_info,
                'fract_preparams' : fract_preparams,
                'fract_full' : fract_full,
                'fract_full_up' : fract_full_up,
                'fract_marker' : fract_marker,
                'fract_marker_up' : fract_marker_up,
                'fract_marker_vis' : fract_marker_vis,
                'fract_mixed_up' : fract_mixed_up,
                'fract_test' : fract_test,
                
                'tp_paths' : tp_paths,
                'tp_tables' : tp_tables,
                'tp_pos' : tp_pos,
                'tp_indata' : tp_indata,
                'tp_data' : tp_data,
                'tp_intermediate' : tp_intermediate,
                'tp_identifiers' : tp_identifiers,
                'tp_info' : tp_info,
                'tp_preparams' : tp_preparams,
                
                'svm_marker' : svm_marker,
                'svm_test' : svm_test,
                'svm_metrics' : svm_metrics,
                
                'marker_sets' : marker_sets,
                'marker_params' : marker_params,
                'marker_conv' : marker_conv,
                'marker_list' : marker_list,
                'marker_fractkey' : marker_fractkey,
                
                'fract_conditions' : fract_conditions,
                
                #'learning_wxyz' : learning_wxyz,
                'learning_xyz' : learning_xyz,
                'results' : results,
                'comparison' : comparison,
                'NN_params' : NN_params,
                
                'status' : status
                
                }
        np.save(filename, file)
    return


def session_open (window, values, filename):
    # filename = sg.popup_get_file('Open Session', no_window=True, file_types=(('Numpy', '*.npy'),))
    # if filename:
    file = np.load(filename, allow_pickle='TRUE').item()
    fract_paths = file['fract_paths']
    fract_tables = file['fract_tables']
    fract_pos = file['fract_pos']
    fract_indata = file ['fract_indata']
    fract_data = file['fract_data']
    fract_std = file['fract_std']
    fract_intermediate = file['fract_intermediate']
    fract_identifiers = file['fract_identifiers']
    fract_info = file['fract_info']
    fract_preparams = file['fract_preparams']
    
    
    
    
    svm_marker = file['svm_marker']
    svm_test = file['svm_test']
    svm_metrics = file['svm_metrics']
    
    fract_full = file['fract_full']
    fract_full_up = file['fract_full_up']
    fract_marker = file['fract_marker']
    fract_marker_up = file['fract_marker_up']
    fract_marker_vis = file['fract_marker_vis']
    fract_mixed_up = file['fract_mixed_up']
    fract_test = file['fract_test']
    
    tp_paths = file['tp_paths']
    tp_tables = file['tp_tables']
    tp_pos = file['tp_pos']
    tp_indata = file['tp_indata']
    tp_data = file['tp_data']
    tp_intermediate = file['tp_intermediate']
    tp_identifiers = file['tp_identifiers']
    tp_info = file['tp_info']
    tp_preparams = file['tp_preparams']
    
    marker_sets = file['marker_sets']
    marker_params = file['marker_params']
    marker_conv = file['marker_conv']
    marker_list = file['marker_list']
    marker_fractkey = file['marker_fractkey']
    
    fract_conditions = file['fract_conditions']
    
    learning_xyz = file['learning_xyz']
    results = file['results']
    comparison = file['comparison']
    NN_params = file['NN_params']
    
    # learning_wxyz = file['learning_wxyz']
    
    status = file['status']
    
    if fract_paths:
        fract_refreshtable(window, fract_tables[fract_paths[0]])
        window['-fractionation_path-'].Update(values = fract_paths, value = fract_paths[0])
    else:
        fract_refreshtable(window, [])
        window['-fractionation_path-'].Update(values = fract_paths, value = '')
    if fract_data['class']:
        fract_buttons(window, True)
    else:
        fract_buttons(window, False)
    if tp_paths:
        tp_refreshtable(window, tp_tables[tp_paths[0]])
        window['-tp_path-'].Update(values = tp_paths, value = tp_paths[0])
    else:
        tp_refreshtable(window, [])
        window['-tp_path-'].Update(values = tp_paths, value = '')
    if tp_data:
        tp_buttons(window, True)
    else:
        tp_buttons(window, False)
    
    if marker_sets:
        refresh_markertable(window, values, marker_sets)
        
        event, values = window.read(timeout = 50)
        #marker_setkey(values, marker_sets)
        #marker_setclass(values, marker_sets)
        refresh_markercols(window, values, marker_sets)
    
    
    
    
        
        
    
    
    
    
    # if marker_list.empty:
    #     CCMPS_actions.enable_markersettings(window, True)
    #     window['-marker_test-'].Update(disabled = False)
    #     window['-marker_profiles-'].Update(disabled = False)
    #     window['-marker_remove-'].Update(disabled = False)
    # else:
    #     CCMPS_actions.enable_markersettings(window, False)
    #     window['-marker_test-'].Update(disabled = True)
    #     window['-marker_profiles-'].Update(disabled = True)
    #     window['-marker_remove-'].Update(disabled = True)
        
    # if fract_data['class'] and not marker_list.empty:
    #     # print('positive')
    #     window['-classification_MOP-'].Update(disabled = False)
    #     #window['-classification_SVM-'].Update(disabled = False)
    # else:
    #     # print('negative')
    #     window['-classification_MOP-'].Update(disabled = True)
        #window['-classification_SVM-'].Update(disabled = True)
        
        #window['-marker_fractkey-'].Update(values = ['[IDENTIFIER]'] + list(fract_info))
    
    # if fract_data['class']:
    #     window['-classification_MOP-'].Update(disabled = False)
    # else:
    #     window['-classification_MOP-'].Update(disabled = True)
    
    # if learning_xyz:
    #     window['-classification_statistics-'].Update(disabled = False)
    # else:
    #     window['-classification_statistics-'].Update(disabled = True)
    
    # if results:
    #     #window['-classification_comparison-'].Update(disabled = False)
    #     window['-status_statistics-'].Update('done!')
    #     window['-export_statistics-'].Update(disabled = False)
    # else:
    #     #window['-classification_comparison-'].Update(disabled = True)
    #     window['-status_statistics-'].Update('missing')
    #     window['-export_statistics-'].Update(disabled = True)
    
    # if comparison:
    #     window['-status_comparison-'].Update('done!')
    #     window['-export_comparison-'].Update(disabled = False)
    # else:
    #     window['-status_comparison-'].Update('missing')
    #     window['-export_comparison-'].Update(disabled = True)
    
    return fract_paths, fract_tables, fract_pos, fract_indata, fract_data, fract_std, fract_intermediate, fract_identifiers, fract_info, fract_preparams, tp_paths, tp_tables, tp_pos, tp_indata, tp_data, tp_intermediate, tp_identifiers, tp_info, tp_preparams, marker_sets, marker_params, marker_conv, fract_conditions, svm_marker, svm_test, svm_metrics, fract_full, fract_full_up, fract_marker, fract_marker_up, fract_marker_vis, fract_mixed_up, fract_test, marker_list, learning_xyz, results, NN_params, comparison, marker_fractkey, status






#------------------------------------------------------------------------------






# def convert_markers(markers, conversion, mode):
#     markerset = pd.DataFrame(columns = ['name'])
#     counter = 1
#     for path in markers:
#         mset = markers[path]['table'][[ markers[path]['identifier_col'] , markers[path]['class_col']         ]]
#         mset.rename(columns = {markers[path]['identifier_col'] : 'name'  , markers[path]['class_col'] : 'class'+str(counter)}, inplace = True)
#         for classname in conversion: 
#             mset['class'+str(counter)].replace({classname : conversion[classname]}, inplace = True)
#             mset['class'+str(counter)].replace(r'^\s*$', np.nan, regex=True, inplace = True)
#             mset = mset[mset['class'+str(counter)].notna()]
#         markerset = pd.merge(markerset, mset, on = 'name', how = 'outer')
#         counter +=1
#     markerset.set_index('name', inplace = True)
#     if mode[0] == 'unite':
#         pass
#     elif mode[0] == 'intersect':
#         markerset.dropna(inplace = True)
#     if mode [1] == 'majority':
#         markerset_final = pd.DataFrame(markerset.mode(axis = 1, dropna = True)[0]).rename(columns = {0 : 'class'})
#     if mode [1] == 'exclude':
#         markerset_final = markerset.mode(axis = 1, dropna = True).fillna(np.nan)
#         if 1 in markerset_final.columns:
#             markerset_final = pd.DataFrame(markerset_final[markerset_final[1].isnull()][0]).rename(columns = {0 : 'class'})
#         else:
#             markerset_final.rename(columns = {0 : 'class'}, inplace = True)
#     return markerset_final



def create_markerprofiles (fract_data, key, fract_info, marker_list):
    profiles = {}
    profiles_vis = {}
    for condition in fract_data['class']:
        profiles[condition] = copy.deepcopy(fract_data['class'][condition])
    for condition in fract_data['vis']:
        profiles_vis[condition] = copy.deepcopy(fract_data['vis'][condition])
    
    fract_marker = {}
    fract_marker_vis = {}
    fract_test = {}
    
    if key == '[IDENTIFIER]':
        for condition in profiles:
            fract_marker[condition] = pd.merge(profiles[condition], marker_list, left_index = True, right_index = True, how = 'left').dropna(subset = ['class'])
            fract_test[condition] = pd.merge(profiles[condition], marker_list, left_index = True, right_index = True, how = 'left')
            fract_test[condition] = fract_test[condition][fract_test[condition]['class'].isna()]
        for condition in profiles_vis:
            fract_marker_vis[condition] = pd.merge(profiles_vis[condition], marker_list, left_index = True, right_index = True)
            
            #profiles_vis[condition] = pd.merge(profiles_vis[condition], marker_list, left_index = True, right_index = True)
            #fract_marker_vis[condition] = pd.merge(profiles_vis[condition], marker_list, left_index = True, right_index = True, how = 'left').dropna(subset = ['class'])
            
    else:
        for condition in profiles:
            profiles[condition] = pd.merge(profiles[condition], fract_info[key].astype(str).applymap(str.upper), left_index=True, right_index=True)
            #profiles[condition] = pd.merge(profiles[condition], fract_info[key].applymap(str.upper), left_index = True, right_index = True)
            
            
            #fract_info_upper = fract_info[key].applymap(str.upper)
            fract_marker[condition] = pd.merge(profiles[condition], marker_list, left_on = key, right_index = True, how = 'left').drop(key, axis = 1).dropna(subset = ['class'])
            fract_test[condition] = pd.merge(profiles[condition], marker_list, left_on = key, right_index = True, how = 'left'). drop(key, axis = 1)
            fract_test[condition] = fract_test[condition][fract_test[condition]['class'].isna()]
        for condition in profiles_vis:
            profiles_vis[condition] = pd.merge(profiles_vis[condition], fract_info[key], left_index = True, right_index = True)
            fract_marker_vis[condition] = pd.merge(profiles_vis[condition], marker_list, left_on = key, right_index = True, how = 'left').drop(key, axis = 1).dropna(subset = ['class'])
        
        
        
        
        # for condition in profiles:
        #     print('check1')
        #     profiles[condition] = pd.merge(profiles[condition], fract_info[key], left_index = True, right_index = True)
        #     print('check2')
        #     profiles_vis[condition] = pd.merge(profiles_vis[condition], fract_info[key], left_index = True, right_index = True)
        #     print('check3')
        #     #print(profiles_vis[condition])
        #     fract_marker[condition] = pd.merge(profiles[condition], marker_list, left_on = key, right_index = True, how = 'left').drop(key, axis = 1).dropna(subset = ['class'])
        #     print('check4')
        #     fract_marker_vis[condition] = pd.merge(profiles_vis[condition], marker_list, left_on = key, right_index = True, how = 'left').drop(key, axis = 1).dropna(subset = ['class'])
        #     print('check5')
        #     fract_test[condition] = pd.merge(profiles[condition], marker_list, left_on = key, right_index = True, how = 'left'). drop(key, axis = 1)
        #     print('check6')
        #     fract_test[condition] = fract_test[condition][fract_test[condition]['class'].isna()]
    
    classnames = {}
    for condition in profiles:
        classnames[condition] = []
        for classname in list(set(fract_marker[condition]['class'].tolist())):
            classnames[condition].append(classname)
            
    
    return fract_marker, fract_marker_vis, fract_test, classnames


# def create_markerprofiles (fract_data, key, fract_info, marker_list):
#     profiles = {}
    
#     if key == '[IDENTIFIER]':
#         for condition in fract_data['class']:
#             profiles[condition] = fract_data['class'][condition]
#     else:
#         for condition in fract_data['class']:
#             profiles[condition] = pd.merge(fract_data['class'][condition], fract_info[key], left_index = True, right_index = True, how = 'left').set_index(key)
    
#     fract_marker = {}
#     fract_test = {}
#     for condition in profiles:
#         fract_marker[condition] = pd.merge(profiles[condition], marker_list, left_index = True, right_index = True, how = 'left').dropna(subset = ['class'])
#         fract_test[condition] = pd.merge(profiles[condition], marker_list, left_index = True, right_index = True, how = 'left')
#         fract_test[condition] = fract_test[condition][fract_test[condition]['class'].isna()]
    
    
#     # for condition in fract_data['class']:
#     #     fract_profiles[condition] = fract_data['class'][condition].set_index(key)
    
    
    
    
#     return fract_marker, fract_test



def upscale (fract_marker, fract_std, key, fract_info, mode):
    stds = {}
    if not key == '[IDENTIFIER]':
        for condition in fract_std['class']:
            stds[condition] = pd.merge(fract_std['class'][condition], fract_info[key], left_index = True, right_index = True, how = 'left').set_index(key)
    
    fract_marker_up = {}
    for condition in fract_marker:
        print('condition', condition)
        
        class_sizes = {}
        for classname in list(set(fract_marker[condition]['class'])):
            class_sizes[classname] = list(fract_marker[condition]['class']).count(classname)
        class_maxsize = max(class_sizes.values())
        
        fract_marker_up[condition] = fract_marker[condition]
        k = 1
        for classname in list(set(fract_marker[condition]['class'])):
            print('class', classname)
            
            data_class_temp = fract_marker[condition].loc[fract_marker[condition]['class'] == classname]
            data_class = data_class_temp.drop(columns = ['class'])
            
            class_difference = abs(class_maxsize - class_sizes[classname])
            # print('maxsize:', class_maxsize)
            # print('class size:', class_sizes[classname])
            # print('difference:', class_difference)
            
            if class_sizes[classname] > class_maxsize:
                ID_rnd = random.sample(list(data_class.index), class_difference -1)
                fract_marker_up[condition].drop(ID_rnd, inplace = True)
            
            if class_sizes[classname] < class_maxsize:
                class_up = pd.DataFrame(columns = data_class.columns)
                
                class_std = data_class.std(axis = 0).to_frame().transpose()
                class_std_flat = class_std.values.flatten()
                # print(class_std)
                
                if mode == 'noised':
                
                
                    for i in range(class_difference):
                        ID_rnd = random.choice(list(data_class.index))
                        name_up = 'up_' + str(k) + '_' + ID_rnd
                        k +=1
                        
                        std_rnd = stds[condition].loc[[ID_rnd]]
                        std_rnd = std_rnd[~std_rnd.index.duplicated(keep = 'first')]
                        
                        profile_rnd = data_class.loc[[ID_rnd]]
                        profile_rnd = profile_rnd[~profile_rnd.index.duplicated(keep = 'first')]
                        
                        
                        
                        list_noised = []
                        
                        for j in range(len(profile_rnd.columns)):
                            col_val = profile_rnd.columns[j]
                            suffix = profile_rnd.columns[j][profile_rnd.columns[j].rfind('_')+1:]
                            col_std = std_rnd.columns[std_rnd.columns.str.endswith(suffix)]
                            sigma = 2* std_rnd[col_std].iloc[0]
                            nv = np.random.normal(profile_rnd[col_val][0], sigma)
                            nv = 0. if nv < 0 else 1. if nv > 1 else nv[0]
                            list_noised.append(nv)
                        
                        profile_noised = pd.DataFrame([list_noised], columns = list(profile_rnd.columns))
                        profile_noised.index = [name_up]
                        profile_noised['class'] = [classname]
                        class_up = class_up.append(profile_noised)
                    fract_marker_up[condition] = fract_marker_up[condition].append(class_up)
                
                if mode == 'average':
                    for i in range(class_difference):
                        ID_rnd_1 = random.choice(list(data_class.index))
                        ID_rnd_2 = random.choice(list(data_class.index))
                        ID_rnd_3 = random.choice(list(data_class.index))
                        name_up = 'up_' + str(k) + '_' + ID_rnd_1 + '_' + ID_rnd_2 + '_' + ID_rnd_3
                        k +=1
                        
                        profile_rnd_1 = data_class.loc[[ID_rnd_1]]
                        profile_rnd_1 = profile_rnd_1[~profile_rnd_1.index.duplicated(keep = 'first')]
                        profile_rnd_2 = data_class.loc[[ID_rnd_2]]
                        profile_rnd_2 = profile_rnd_2[~profile_rnd_2.index.duplicated(keep = 'first')]
                        profile_rnd_3 = data_class.loc[[ID_rnd_3]]
                        profile_rnd_3 = profile_rnd_3[~profile_rnd_3.index.duplicated(keep = 'first')]
                        
                        profile_av = pd.concat([profile_rnd_1, profile_rnd_2, profile_rnd_3]).median(axis = 0).to_frame().transpose()
                        
                        profile_av.index = [name_up]
                        profile_av['class'] = [classname]
                        class_up = class_up.append(profile_av)
                    fract_marker_up[condition] = fract_marker_up[condition].append(class_up)
                
                if mode == 'noisedaverage':
                    for i in range(class_difference):
                        # print(i)
                        ID_rnd_1 = random.choice(list(data_class.index))
                        ID_rnd_2 = random.choice(list(data_class.index))
                        ID_rnd_3 = random.choice(list(data_class.index))
                        name_up = 'up_' + str(k) + '_' + ID_rnd_1 + '_' + ID_rnd_2 + '_' + ID_rnd_3
                        k +=1
                        
                        # class_std = data_class.std(axis = 1).to_frame().transpose()
                        
                        profile_rnd_1 = data_class.loc[[ID_rnd_1]]
                        profile_rnd_1 = profile_rnd_1[~profile_rnd_1.index.duplicated(keep = 'first')]
                        profile_rnd_2 = data_class.loc[[ID_rnd_2]]
                        profile_rnd_2 = profile_rnd_2[~profile_rnd_2.index.duplicated(keep = 'first')]
                        profile_rnd_3 = data_class.loc[[ID_rnd_3]]
                        profile_rnd_3 = profile_rnd_3[~profile_rnd_3.index.duplicated(keep = 'first')]
                        
                        profile_av = pd.concat([profile_rnd_1, profile_rnd_2, profile_rnd_3]).median(axis = 0).to_frame().transpose()
                        # print(len(profile_av))
                        profile_av_flat = profile_av.values.flatten()
                        
                        list_noised = []
                        
                        # print(class_std)
                        # 
                        # for j in range(len(class_std.columns)):
                        #     # sigma = 2*class_std[class_std.columns[j]]
                        #     nv = np.random.normal(profile_av_flat, 2* class_std_flat, size = profile_av.shape)
                        #     nv = 0. if nv < 0 else 1. if nv > 1 else nv[0]
                        #     # print(sigma)
                        nv = np.random.normal(profile_av_flat, 2* class_std_flat, size = profile_av.shape)
                        nv = np.where(nv > 1, 1, np.where(nv < 0, 0, nv))
                        
                        # print(nv)
                        
                        # values = np.where(values > 1, 1, np.where(values < 0, 0, values))
                        
                        
                        
                        profile_noised = pd.DataFrame(nv, columns = profile_av.columns)
                        profile_noised.index = [name_up]
                        profile_noised['class'] = [classname]
                        class_up = class_up.append(profile_noised)
                        
                        # profile_av.index = [name_up]
                        # profile_av['class'] = [classname]
                        # class_up = class_up.append(profile_av)
                    fract_marker_up[condition] = fract_marker_up[condition].append(class_up)
                        
                        
                    
        #noised_df = pd.DataFrame([noised_values], columns=df_mean.columns)
        
        
        
        
        
        
    return fract_marker_up, class_maxsize






def marker_mix (fract_marker_up):
    fract_mixed_up = {}
    
    for condition in fract_marker_up:
        class_list = list(set(list(fract_marker_up[condition]['class'])))
        combinations = [(a,b) for idx, a in enumerate(class_list) for b in class_list[idx + 1:]]
        
        
        fract_mixed_up[condition] = pd.DataFrame(columns = fract_marker_up[condition].drop(columns = ['class']).columns)
        for classname in class_list:
            fract_mixed_up[condition][classname] = 0.
        
        
        
        
        cur = 1
        for comb in combinations:
            
            
            profiles_own = fract_marker_up[condition].copy().loc[fract_marker_up[condition]['class'] == comb[0]].drop(columns = ['class'])
            profiles_other = fract_marker_up[condition].copy().loc[fract_marker_up[condition]['class'] == comb[1]].drop(columns = ['class'])
            
            new_index = [f"{i}_{j}" for i, j in zip(profiles_own.index, profiles_other.index)]
            mix_steps = [i/10 for i in range(0, 11)]
            
            for part in mix_steps:
                
                new_index_part = [f"{i + cur}_{value}" for i, value in enumerate(new_index)]
                own_part = profiles_own.multiply(part)
                other_part = profiles_other.multiply((1 - part))
                
                own_part.index = new_index_part
                other_part.index = new_index_part
                
                #profiles_mixed = profiles_own.add(profiles_other, fill_value = 0)
                profiles_mixed = own_part + other_part
                
                for classname in class_list:
                    if classname == comb[0]:
                        profiles_mixed[classname] = part
                    elif classname == comb[1]:
                        profiles_mixed[classname] = (1 - part)
                    else:
                        profiles_mixed[classname] = 0.
                
                
                
                fract_mixed_up[condition] = pd.concat([fract_mixed_up[condition], profiles_mixed])
                
                cur += len(profiles_mixed)
        
        # fract_marker_dummies = pd.get_dummies(fract_marker_up['Class'])
        # fract_mixed_up = pd.concat([])
        
        # df = pd.concat([df, class_dummies], axis=1).drop(columns=['Class'])

        
        #fract_mixed_up[condition] = pd.concat([fract_mixed_up[condition], fract_marker_up[condition]])
    return fract_mixed_up





# res = [(a, b) for idx, a in enumerate(test_list) for b in test_list[idx + 1:]]











# df_add = df1.add(df2, fill_value=0)

# df.loc[df['column_name'] == some_value]

# def marker_mix (fract_marker, class_maxsize):
#     mix_steps = [i/5 for i in range(1, 5)]
#     mix_steps_reverse = [1-i for i in mix_steps]
    
#     fract_marker_mixed = {}
    
#     for condition in fract_marker:
#         # print('mixing classes for', condition)
        
#         fract_marker_mixed[condition] = pd.DataFrame(columns = fract_marker[condition].columns)
#         classlist = list(set(fract_marker[condition]['class']))
#         for newclass in classlist:
#             fract_marker_mixed[condition][newclass] = 0
        
#         k = 1
#         for class_own in list(set(fract_marker[condition]['class'])):
#             classlist.remove(class_own)
#             data_class_temp_own = fract_marker[condition].loc[fract_marker[condition]['class'] == class_own]
#             data_class_own = data_class_temp_own.drop(columns = ['class'])
            
#             for class_other in classlist:
#                 data_class_temp_other = fract_marker[condition].loc[fract_marker[condition]['class'] == class_other]
#                 data_class_other = data_class_temp_other.drop(columns = ['class'])
                
#                 print('mixing classes for', condition, class_own, class_other)
                
#                 for i in range(int(0.1*class_maxsize)):
#                     ID_rnd_own = random.choice(list(data_class_own.index))
#                     ID_rnd_other = random.choice(list(data_class_other.index))
#                     # print('own ID:', ID_rnd_own)
#                     # print('other ID:', ID_rnd_other)
                    
#                     profile_rnd_own = data_class_own.loc[[ID_rnd_own]]
#                     profile_rnd_own = profile_rnd_own[~profile_rnd_own.index.duplicated(keep = 'first')]
#                     profile_rnd_other = data_class_other.loc[[ID_rnd_other]]
#                     profile_rnd_other = profile_rnd_other[~profile_rnd_other.index.duplicated(keep = 'first')]
                    
#                     for r in range(len(mix_steps)):
#                         ID_mixed = 'mixed_' + str(k) + '_' + str(mix_steps[r]) + '*' + ID_rnd_own + '_' + str(mix_steps_reverse[r]) + '*' + ID_rnd_other
                        
#                         profile_rnd_own *= mix_steps[r]
#                         profile_rnd_other *= mix_steps_reverse[r]
                        
#                         profile_rnd_own.index = [ID_mixed]
#                         profile_rnd_other.index = [ID_mixed]
                        
#                         profile_rnd_combined = profile_rnd_own.add(profile_rnd_other, fill_value = 0)
#                         # profile_rnd_combined = profile_rnd_combined.div(profile_rnd_combined.sum(axis = 1), axis = 0)
                        
#                         for newclass in list(set(fract_marker[condition]['class'])):
#                             profile_rnd_combined[newclass] = 0.
                        
#                         # print(profile_rnd_combined)
#                         profile_rnd_combined[class_own] = [mix_steps[r]]
#                         profile_rnd_combined[class_other] = [mix_steps_reverse[r]]
#                         profile_rnd_combined.index = [ID_mixed]
                        
#                         fract_marker_mixed[condition] = pd.merge(fract_marker_mixed[condition], profile_rnd_combined, left_index = True, right_index = True, how = 'outer')
                        
#                         k +=1
                        
                        
#     return fract_marker_mixed











# concatenated_df = pd.concat([df1, df2, df3])

# # Calculate the median for each column
# median_values = concatenated_df.median()

# # Create a new dataframe with the median values
# median_df = pd.DataFrame(median_values).transpose()








def create_fullprofiles (fract_marker, fract_test, fract_marker_up):
    fract_full = {}
    fract_full_up = {}
    for condition in fract_test:
        fract_full[condition] = pd.concat([fract_test[condition], fract_marker[condition]])
        fract_full_up[condition] = pd.concat([fract_test[condition], fract_marker_up[condition]])
    
    
    # fract_full = pd.concat([fract_test, fract_marker])
    # fract_full_up = pd.concat([fract_test, fract_marker_up])
    
    
    return fract_full, fract_full_up





# suffix = sample[sample.rfind('_')+1:]

# df.columns[df.columns.str.endswith("_Emp")]



def enable_markersettings (window, is_enabled):
    for element in ['-marker_add-', '-marker_remove-', '-marker_key-', '-marker_class-', '-marker_fractkey-', '-marker_parameters-', '-marker_manage-', '-marker_accept-', '-marker_profiles-', '-marker_test-']:
        window[element].Update(disabled = not is_enabled)
    
    
    for element in ['-marker_reset-']:
        window[element].Update(disabled = is_enabled)
    
    
    # window['-marker_reset-'].Update(disabled = is_enabled)
    if is_enabled:
        window['-status_marker-'].Update(value = 'missing', text_color = 'white')
    else:
        window['-status_marker-'].Update(value = 'ready!', text_color = 'dark green')
    return



# data_up = data_learning
# n = 1
# for organelle in list(set(data_learning['class'])):
#     data_class_temp = data_learning.loc[data_learning['class'] == organelle]
#     data_class = data_class_temp.drop(columns=['class'])
    
    
#     #data_class = data_learning.loc[data_learning['class'] == organelle].drop(columns=['class'])
    
    
#     class_difference = abs(batch_target - class_sizes[organelle])
#     print(class_difference)
    
#     if class_sizes[organelle] > batch_target:
#         ID_rnd = random.sample(list(data_class.index), class_difference-1)
#         data_up.drop(ID_rnd, inplace = True)
    
#     if class_sizes[organelle] < batch_target:
#         class_up = pd.DataFrame(columns = data_class.columns)
#         for i in range(class_difference):
#             ID_rnd = random.choice(list(data_class.index))
#             name_up = 'up_' + str(n) + '_' + ID_rnd
#             n = n +1
            
#             std_rnd = data_std.loc[[ID_rnd]]
            
#             profile_rnd = data_class.loc[[ID_rnd]]
#             profile_rnd = profile_rnd[~profile_rnd.index.duplicated(keep='first')]
            
#             list_noised = []
            
#             for j in range(len(profile_rnd.columns)):
#                 col_val = profile_rnd.columns[j]
#                 col_std = std_rnd.columns[j]
#                 sigma = 0.5* std_rnd[col_std][0]
#                 nv = np.random.normal(profile_rnd[col_val][0], sigma)
#                 nv = 0. if nv < 0 else 1. if nv > 1 else nv
#                 list_noised.append(nv)
            
#             profile_noised = pd.DataFrame([list_noised], columns= list(profile_rnd.columns))
#             profile_noised.index = [name_up]
#             profile_noised['class'] = [organelle]
#             class_up = class_up.append(profile_noised)
        
#         data_up = data_up.append(class_up)

# class_sizes = {}
# for organelle in list(set(data_learning['class'])):
#     class_sizes[organelle] = list(data_learning['class']).count(organelle)
# class_maxsize = max(class_sizes.values())













def default_status ():
    status = {'fractionation_data' : False,
              'tp_data' : False,
              'lipidome_data' : False,
              'lipidome_total' : False,
              'marker_file' : False,
              'marker_matched' : False,
              'training' : False,
              'proteome_prediction' : False,
              'lipidome_prediction' : False,
              'comparison_global' : False,
              'comparison_class' : False}
    return status









def refresh_window (window, status):
    for element in ['-fractionation_reset-', '-fractionation_summary-']:
        window[element].Update(disabled = not status['fractionation_data'])
    for element in ['-fractionation_add-', '-fractionation_remove-',
                    '-fractionation_edit_remove-', '-fractionation_edit_keep-', '-fractionation_edit_condition-',
                    '-fractionation_edit_replicate-', '-fractionation_edit_fractions-', '-fractionation_edit_identifier-',
                    '-fractionation_parameters-', '-fractionation_start-']:
        window[element].Update(disabled = status['fractionation_data'])
    
    for element in ['-tp_reset-', '-tp_summary-', '-tp_export-']:
        window[element].Update(disabled = not status['tp_data'])
    for element in ['-tp_add-', '-tp_remove-',
                    '-tp_edit_remove-', '-tp_edit_keep-', '-tp_edit_condition-', '-tp_edit_identifier-',
                    '-tp_parameters-', '-tp_start-']:
        window[element].Update(disabled = status['tp_data'])
    
    
    for element in ['-marker_remove-', '-marker_manage-', '-marker_accept-']:
        if status['marker_matched']:
            window[element].Update(disabled = True)
        else:
            window[element].Update(disabled = not status['marker_file'])
    for element in ['-marker_test-', '-marker_profiles-']:
        window[element].Update(disabled = not status['marker_file'])
    for element in ['-marker_reset-']:
        window[element].Update(disabled = not status['marker_matched'])
    for element in ['-marker_add-', '-marker_parameters-', '-marker_preset-']:
        window[element].Update(disabled = status['marker_matched'])
    
    for element in ['-statistic_import-']:
        window[element].Update(disabled = status['comparison_global'])
    
    if status['fractionation_data']:
        window['-status_fract-'].Update('ready', text_color = 'dark green')
    else:
        window['-status_fract-'].Update('none', text_color = 'black')
    if status['tp_data']:
        window['-status_tp-'].Update('ready', text_color = 'dark green')
    else:
        window['-status_tp-'].Update('none', text_color = 'black')
    if status['lipidome_data']:
        window['-status_fract_lipid-'].Update('ready', text_color = 'dark green')
    else:
        window['-status_fract_lipid-'].Update('none', text_color = 'black')
    if status['marker_matched']:
        window['-status_marker-'].Update('ready', text_color = 'dark green')
    else:
        window['-status_marker-'].Update('none', text_color = 'black')
    if status['lipidome_total']:
        window['-status_total_lipid-'].Update('ready', text_color = 'dark green')
    else:
        window['-status_total_lipid-'].Update('none', text_color = 'black')
    
    
    for element in ['-classification_MOP-']:
        if status['fractionation_data'] and status['marker_matched']:
            window[element].Update(disabled = status['training'])
        else:
            window[element].Update(disabled = True)
    
    for element in ['-classification_validation-', '-classification_reset-']:
        window[element].Update(disabled = not status['training'])
    
    
    
    if status['training']:
        for element in ['-statistic_predict-']:
            window[element].Update(disabled = status['proteome_prediction'])
        for element in ['-statistic_export-', '-statistic_report-',
                        '-statistic_reset-', '-statistic_heatmap-', '-statistic_distribution-',]:
            window[element].Update(disabled = not status['proteome_prediction'])
        
        if status['proteome_prediction']:
            for element in ['-global_run-']:
                window[element].Update(disabled = status['comparison_global'])
            for element in ['-global_heatmap-', '-global_distance-',
                            '-global_report-', '-global_reset-']:
                window[element].Update(disabled = not status['comparison_global'])
            
            if status['comparison_global'] and status['tp_data']:
                for element in ['-class_run-']:
                    window[element].Update(disabled = status['comparison_class'])
                for element in ['-class_heatmap-', '-class_reorganization-',
                                '-class_report-', '-class_reset-']:
                    window[element].Update(disabled = not status['comparison_class'])
            else:
                for element in ['-class_run-',
                                '-class_heatmap-', '-class_reorganization-',
                                '-class_report-', '-class_reset-']:
                    window[element].Update(disabled = True)
                
                if status['lipidome_data']:
                    for element in ['-lipidome_predict-']:
                        window[element].Update(disabled = status['lipidome_prediction'])
                    for element in ['-lipidome_report-', '-lipidome_reset-',
                                    '-lipidome_heatmap-', '-lipidome_reorganization-',
                                    '-lipidome_density-', '-lipidome_composition-']:
                        window[element].Update(disabled = not status['lipidome_prediction'])
                else:
                    for element in ['-lipidome_predict-',
                                    '-lipidome_report-', '-lipidome_reset-',
                                    '-lipidome_heatmap-', '-lipidome_reorganization-',
                                    '-lipidome_density-', '-lipidome_composition-']:
                        window[element].Update(disabled = True)
        
        else:
            for element in ['-lipidome_predict-',
                            '-lipidome_report-', '-lipidome_reset-',
                            '-lipidome_heatmap-', '-lipidome_reorganization-',
                            '-lipidome_density-', '-lipidome_composition-',
                            
                            '-global_run-', '-global_heatmap-', '-global_distance-',
                            '-global_report-', '-global_reset-',
                            
                            '-class_run-', '-class_heatmap-', '-class_reorganization-',
                            '-class_report-', '-class_reset-']:
                window[element].Update(disabled = True)
        
    else:
        for element in ['-statistic_predict-',
                        '-statistic_export-', '-statistic_report-',
                        '-statistic_reset-', '-statistic_heatmap-', '-statistic_distribution-',
                        
                        '-lipidome_predict-',
                        '-lipidome_report-', '-lipidome_reset-', '-lipidome_heatmap-',
                        '-lipidome_reorganization-', '-lipidome_density-', '-lipidome_composition-',
                        
                        '-global_run-', '-global_heatmap-', '-global_distance-',
                        '-global_report-', '-global_reset-',
                        
                        '-class_run-', '-class_heatmap-', '-class_reorganization-',
                        '-class_report-', '-class_reset-']:
            window[element].Update(disabled = True)
    
    
    
    
    
    
    
    return























