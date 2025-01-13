## SHOW MARKERS



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

def create_line_plot(data, title=None):
    fig, ax = plt.subplots(figsize=(8, 6))  # Adjust the figure size as needed
    for column in data.columns:
        ax.plot(data.index, data[column], label=column)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)  # Place the legend below the plot
    ax.set_xlabel('fractions')
    ax.set_ylabel('normalized intensity')
    ax.set_xticks([])  # Remove x-tick labels
    ax.set_yticks([0, 1])
    ax.set_yticklabels(['0', '1'])
    ax.set_xlim(0, len(data.index) - 1)  # Set x-axis limits
    if title:
        plt.title(title)
    plt.ylim(0, 1)
    fig.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust layout to make room for the legend
    return fig

def update_class_info(marker_list, classnames, data):
    class_info = []
    for classname in classnames:
        count = data[data.index.isin(marker_list[marker_list['class'] == classname].index)].shape[0]
        class_info.append([classname, count])
    return class_info

def SM_exec(fract_data, fract_info, marker_list, key):
    profiles_dict = {}
    class_info_dict = {}
    distinct_profiles_dict = {}
    
    for condition in fract_data['vis']:
        data = pd.merge(fract_data['vis'][condition], fract_info[key], left_index=True, right_index=True, how='left')
        data.set_index(key, inplace=True)
        
        median_classprofiles = {}
        distinct_profiles = {}
        
        classnames = list(set(marker_list['class']))
        for classname in classnames:
            marker_class = marker_list[marker_list['class'] == classname]
            data_class = data[data.index.isin(marker_class.index)]
            median_classprofiles[classname] = data_class.median()
            distinct_profiles[classname] = data_class
        
        profiles_df = pd.DataFrame(median_classprofiles)
        profiles_dict[condition] = profiles_df
        class_info_dict[condition] = update_class_info(marker_list, classnames, data)
        distinct_profiles_dict[condition] = distinct_profiles
    
    sg.theme('DarkTeal11')
    
    condition = list(profiles_dict.keys())[0]
    layout_SM = [
        [sg.Text('Select Condition:'), sg.Combo(list(profiles_dict.keys()), key='-condition-', enable_events=True, default_value=condition, readonly=True, size = (25,1))],
        [sg.Canvas(key='-CANVAS-'), sg.Table(values=class_info_dict[condition], headings=['Class', 'n'], key='-CLASSINFO-', col_widths=[20, 10], auto_size_columns=False, justification='left', num_rows=26)],
        [sg.Button('Export all Conditions...', key='-EXPORT-', size = (20,1))]
    ]
    
    window_SM = sg.Window('Marker Parameters', layout_SM, finalize=True, size=(920, 520))  # Adjust window width
    
    # Initial drawing
    fig = create_line_plot(profiles_dict[condition], title=condition)
    figure_agg = draw_figure(window_SM['-CANVAS-'].TKCanvas, fig)
    
    while True:
        event_SM, values_SM = window_SM.read()
        
        if event_SM == sg.WIN_CLOSED:
            break
        elif event_SM == '-condition-':
            condition = values_SM['-condition-']
            fig = create_line_plot(profiles_dict[condition], title=condition)
            figure_agg = update_figure(window_SM['-CANVAS-'].TKCanvas, figure_agg, fig)
            window_SM['-CLASSINFO-'].update(values=class_info_dict[condition])
        elif event_SM == '-EXPORT-':
            folder_path = sg.popup_get_folder('Select Folder')
            if folder_path:
                # Save the main Excel file with all conditions and median profiles
                with pd.ExcelWriter(os.path.join(folder_path, 'markerprofiles_combined.xlsx')) as writer:
                    for cond, df in profiles_dict.items():
                        df.to_excel(writer, sheet_name=cond)
                
                # Save individual Excel files for each condition with distinct profiles for each class
                for cond, distinct_profiles in distinct_profiles_dict.items():
                    with pd.ExcelWriter(os.path.join(folder_path, f'markerprofiles_{cond}.xlsx')) as writer:
                        for classname, df in distinct_profiles.items():
                            df.to_excel(writer, sheet_name=classname)
                
                # Save the plot
                for cond, df in profiles_dict.items():
                    fig = create_line_plot(df, title=cond)
                    fig.savefig(os.path.join(folder_path, f'markerprofiles_{cond}.pdf'), format='pdf')
                    plt.close(fig)
    
    window_SM.close()
    return

# Example usage
# fract_data, fract_info, marker_list, key = your data here
# profiles_dict = SM_exec(fract_data, fract_info, marker_list, key)


























