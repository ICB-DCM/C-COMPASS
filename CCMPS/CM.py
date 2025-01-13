##### CLASS MANAGER

import PySimpleGUI as sg
import numpy as np
import pandas as pd






def refresh_conversion(conversion, values):
    for o in conversion:
        if values['--'+o+'--'] == True:
            conversion[o] = values['--'+o+'_class--']
        elif values['--'+o+'--'] == False:
            conversion[o] = np.nan
    return conversion
    
def CM_exec (markers, marker_conv):
    sg.theme('DarkTeal11')
    conv_old = marker_conv
    
    num_names = 0
    for name in marker_conv:
        if not pd.isnull(marker_conv[name]):
            num_names += 1
    
    layout_column = [
        [sg.Text('Annotation', size = (25,1)),
         sg.Text('Class', size = (20,1))],
        [sg.Text('--------------------------------------------------------------------------------')],
        *[[sg.Checkbox(o, default = not pd.isnull(marker_conv[o]), enable_events = True, size = (20,5), key = '--'+o+'--'),
           sg.InputText(str(marker_conv[o]), visible = not pd.isnull(marker_conv[o]), size = (20,5), key = '--'+o+'_class--')] for o in marker_conv]
        ]
    
    layout_CM = [
        [sg.Frame(layout = [
                      [sg.Column(layout = layout_column, size = (380,340), scrollable = True, vertical_scroll_only = True)
                      ]
                ],
                title = 'Classes', size = (400,380)),
        
        sg.Column(layout = [
                      [sg.Text('initial annotations: '),
                       sg.Text(str(len(marker_conv)))],
                      [sg.Text('used annotations: '),
                       sg.Text(str(num_names), key = '-num_anno-')],
                      [sg.Button('Accept', size = (15,1), enable_events = True, button_color = 'dark green', key = '--accept--')],
                      [sg.Button('Cancel', size = (15,1), enable_events = True, button_color = 'black', key = '--cancel--')]
                ], size = (200,340))]
        ]
    
    window_CM = sg.Window('Classes', layout_CM, size = (600,400))
    
    while True:
        event_CM, values_CM = window_CM.read()
        for k in marker_conv:
            if event_CM == '--'+k+'--':
                window_CM['--'+k+'_class--'].Update(visible = values_CM['--'+k+'--'])
                if values_CM['--'+k+'--'] == True:
                    window_CM['--'+k+'_class--'].Update(value = k)
                    num_names += 1
                    window_CM['-num_anno-'].Update (value = str(num_names))
                elif values_CM['--'+k+'--'] == False:
                    window_CM['--'+k+'_class--'].Update(value = False)
                    num_names -= 1
                    window_CM['-num_anno-'].Update (value = str(num_names))
        
        if event_CM == sg.WIN_CLOSED or event_CM == '--cancel--':
            marker_conv = conv_old
            break
        if event_CM == '--accept--':
            marker_conv = refresh_conversion(marker_conv, values_CM)
            break
        
    window_CM.close()
    return marker_conv


