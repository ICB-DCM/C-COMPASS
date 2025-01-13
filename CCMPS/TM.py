##### TEST MARKER

import PySimpleGUI as sg
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os

def draw_figure(canvas, figure):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg

def update_figure(canvas, figure_agg, figure):
    figure_agg.get_tk_widget().forget()
    plt.close('all')
    return draw_figure(canvas, figure)

def create_heatmap(dataframe, title=None):
    fig, ax = plt.subplots(figsize=(8, 8))  # Adjust the figure size as needed
    cax = ax.matshow(dataframe, cmap='coolwarm', vmin=-1, vmax=1)
    fig.colorbar(cax)
    ax.set_xticks(range(len(dataframe.columns)))
    ax.set_xticklabels(dataframe.columns, rotation=90, fontsize=8)  # Rotate x-axis labels 90 degrees
    ax.set_yticks(range(len(dataframe.index)))
    ax.set_yticklabels(dataframe.index, fontsize=8)
    plt.subplots_adjust(top=0.8, bottom=0.1, left=0.2)  # Adjust the top and bottom margins
    if title:
        plt.title(title)
    return fig

def update_class_info(marker_list, classnames, data):
    class_info = []
    for classname in classnames:
        count = data[data.index.isin(marker_list[marker_list['class'] == classname].index)].shape[0]
        class_info.append([classname, count])
    return class_info

def TM_exec(fract_data, fract_info, marker_list, key):
    correlation_matrices = {}
    class_info_dict = {}
    
    for condition in fract_data['vis']:
        data = pd.merge(fract_data['vis'][condition], fract_info[key], left_index=True, right_index=True, how='left')
        data.set_index(key, inplace=True)
        
        median_classprofiles = {}
        
        classnames = list(set(marker_list['class']))
        for classname in classnames:
            marker_class = marker_list[marker_list['class'] == classname]
            data_class = data[data.index.isin(marker_class.index)]
            median_classprofiles[classname] = data_class.median().to_numpy()
        
        correlation_matrix = np.zeros((len(classnames), len(classnames)))
        
        for i, class1 in enumerate(classnames):
            for j, class2 in enumerate(classnames):
                if i == j:
                    correlation_matrix[i, j] = 1.0
                else:
                    correlation_matrix[i, j] = np.corrcoef(median_classprofiles[class1], median_classprofiles[class2])[0, 1]
        
        correlation_df = pd.DataFrame(correlation_matrix, index=classnames, columns=classnames)
        correlation_matrices[condition] = correlation_df
        class_info_dict[condition] = update_class_info(marker_list, classnames, data)
    
    sg.theme('DarkTeal11')
    
    condition = list(correlation_matrices.keys())[0]
    layout_TM = [
        [sg.Text('Select Condition:'), sg.Combo(list(correlation_matrices.keys()), key='-condition-', enable_events=True, default_value=condition, readonly=True, size = (25,1))],
        [sg.Canvas(key='-CANVAS-'), sg.Table(values=class_info_dict[condition], headings=['Class', 'n'], key='-CLASSINFO-', col_widths=[20, 10], auto_size_columns=False, justification='left', num_rows=35)],
        [sg.Button('Export all Conditions...', key='-EXPORT-', size = (20,1))]
    ]
    
    window_TM = sg.Window('Marker Parameters', layout_TM, finalize=True, size=(900, 650))  # Increase window width
    
    # Initial drawing
    fig = create_heatmap(correlation_matrices[condition], title=condition)
    figure_agg = draw_figure(window_TM['-CANVAS-'].TKCanvas, fig)
    
    while True:
        event_TM, values_TM = window_TM.read()
        
        if event_TM == sg.WIN_CLOSED:
            break
        elif event_TM == '-condition-':
            condition = values_TM['-condition-']
            fig = create_heatmap(correlation_matrices[condition], title=condition)
            figure_agg = update_figure(window_TM['-CANVAS-'].TKCanvas, figure_agg, fig)
            window_TM['-CLASSINFO-'].update(values=class_info_dict[condition])
        elif event_TM == '-EXPORT-':
            folder_path = sg.popup_get_folder('Select Folder')
            if folder_path:
                for cond, df in correlation_matrices.items():
                    # Save the plot
                    fig = create_heatmap(df, title=cond)
                    fig.savefig(os.path.join(folder_path, f'{cond}.pdf'), format='pdf')
                    plt.close(fig)
                
                # Save all data to an Excel file
                with pd.ExcelWriter(os.path.join(folder_path, 'correlation_matrices.xlsx')) as writer:
                    for cond, df in correlation_matrices.items():
                        df.to_excel(writer, sheet_name=cond)
    
    window_TM.close()
    return

# Example usage
# fract_data, fract_info, marker_list, key = your data here
# correlation_matrices = TM_exec(fract_data, fract_info, marker_list, key)











































