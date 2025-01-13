##### PREPARAMETERS

import PySimpleGUI as sg
import PPMS_actions as action
import copy

def fract_default ():
    params_default = {
        'class' : {
            'scale1' : [
                True,
                'area',
                ],
            'corrfilter' : False,
            'scale2' : [
                False,
                'area'
                ],
            'zeros' : True,
            'combination' : 'separate'
            },
        'vis' : {
            'scale1' : [
                True,
                'minmax',
                ],
            'corrfilter' : False,
            'scale2' : [
                True,
                'minmax'
                ],
            'zeros' : True,
            'combination': 'median'
            },
        'global' : {
            'missing' : [
                True,
                '1'
                ],
            'minrep' : [
                True,
                '2'
                ],
            'outcorr' : False
            }
        }
    return params_default

def tp_default ():
    params_default = {
        'minrep' : 2,
        'imputation' : 'normal'
        }
    return params_default


def NN_default ():
    params_default = {
        'upsampling' : True,
        'upsampling_method' : 'noisedaverage',
        'upsampling_noise' : 2,
        'AE' : 'none',
        'AE_activation' : 'leakyrelu',
        'AE_out' : 'sigmoid',
        'AE_epochs' : 20,
        'svm_filter' : False,
        'mixed_part' : 4,
        'mixed_batch' : 0.05,
        'NN_optimization' : 'long',
        'NN_activation' : 'relu',
        'class_activation' : 'linear',
        'class_loss' : 'mean_squared_error',
        'regularization' : 'none',
        'optimizers' : ['adam', 'rmsprop', 'sgd'],
        'NN_epochs' : 20,
        'rounds' : 3,
        'subrounds' : 10,
        'reliability' : 95
        }
    return params_default



def PPMS_exec_fract (params_old):
    combination_both = action.set_combinations(params_old)
    
    sg.theme('DarkTeal11')
    
    layout_class = [[sg.Checkbox('pre-scaling', key = '--class_scale1--', disabled = False, enable_events = True, default = params_old['class']['scale1'][0]),
                     sg.Combo(['minmax', 'area'], key = '--class_scale1_how--', size = (10,1), disabled = not params_old['class']['scale1'][0], enable_events = False, readonly = True, default_value = params_old['class']['scale1'][1])
                    ],
                    [sg.Checkbox('exclude proteins from worst correlated replicate', key = '--class_corrfilter--', disabled = False, enable_events = False, default = params_old['class']['corrfilter'])],
                    [sg.Checkbox('median profile', key = '--class_median--', disabled = False, enable_events = True, default = combination_both['class_median'], tooltip = ' for consistent/reproducible replicates '),
                     sg.Checkbox('concatenated profiles', key = '--class_concat--', disabled = False, enable_events = True, default = combination_both['class_concat'], tooltip = ' for variations over replicate '),
                     sg.Checkbox('process separately', key = '--class_separate--', disabled = False, enable_events = True, default = combination_both['class_separate'], tooltip = ' for other purposes ')],
                    [sg.Checkbox('post-scaling', key = '--class_scale2--', disabled = False, enable_events = True, default = params_old['class']['scale2'][0]),
                     sg.Combo(['minmax', 'area'], key = '--class_scale2_how--', size = (10,1), disabled = not params_old['class']['scale2'][1], enable_events = False, readonly = True, default_value = params_old['class']['scale2'][1])],
                    [sg.Checkbox('remove baseline profiles (zeroes)', key = '--class_zeros--', disabled = False, enable_events = False, default = params_old['class']['zeros'])],
                    ]
    layout_vis = [[sg.Checkbox('pre-scaling', key = '--vis_scale1--', disabled = False, enable_events = True, default = params_old['vis']['scale1'][0]),
                     sg.Combo(['minmax', 'area'], key = '--vis_scale1_how--', size = (10,1), disabled = not params_old['vis']['scale1'][0], enable_events = False, readonly = True, default_value = params_old['vis']['scale1'][1])
                    ],
                    [sg.Checkbox('exclude proteins from worst correlated replicate', key = '--vis_corrfilter--', disabled = False, enable_events = False, default = params_old['vis']['corrfilter'])],
                    [sg.Checkbox('median profile', key = '--vis_median--', disabled = False, enable_events = True, default = combination_both['vis_median'], tooltip = ' for consistent/reproducible replicates '),
                     sg.Checkbox('concatenated profiles', key = '--vis_concat--', disabled = False, enable_events = True, default = combination_both['vis_concat'], tooltip = ' for variations over replicate '),
                     sg.Checkbox('process separately', key = '--vis_separate--', disabled = False, enable_events = True, default = combination_both['vis_separate'], tooltip = ' for other purposes ')],
                    [sg.Checkbox('post-scaling', key = '--vis_scale2--', disabled = False, enable_events = True, default = params_old['vis']['scale2'][0]),
                     sg.Combo(['minmax', 'area'], key = '--vis_scale2_how--', size = (10,1), disabled = not params_old['vis']['scale2'][0], enable_events = False, readonly = True, default_value = params_old['vis']['scale2'][1])],
                    [sg.Checkbox('remove baseline profiles (zeroes)', key = '--vis_zeros--', disabled = False, enable_events = False, default = params_old['class']['zeros'])],
                    ]
    
    layout_PPMS = [
                [sg.TabGroup([[sg.Tab(' - Classification - ', layout_class), sg.Tab(' - Visualization - ', layout_vis)]],
                tab_location = 'top', tab_background_color = 'grey', size = (450,225)),
                sg.Frame(layout = [
                    [sg.Checkbox('min. valid fractions:  ', key = '--glob_missing--', disabled = False, enable_events = True, default = params_old['global']['missing'][0])],
                    [sg.Spin(values=(list(range(1,100))), size = (10,2), key = '--glob_minval--', disabled = not params_old['global']['missing'][0], enable_events = False, readonly = True, initial_value = params_old['global']['missing'][1], text_color = 'black')],
                    [sg.Text('-' *45)],
                    [sg.Checkbox('found in at least: ', key = '--glob_minrep--', disabled = False, enable_events = True, default = params_old['global']['minrep'][0])],
                    [sg.Spin(values=(list(range(1,10))), size = (10,2), key = '--glob_mincount--', disabled = not params_old['global']['minrep'][0], enable_events = False, readonly = True, initial_value = params_old['global']['minrep'][1], text_color = 'black'),
                     sg.Text(' replicates', text_color = 'light grey')],
                    [sg.Text('-' *45)],
                    [sg.Checkbox('calculate outer correlations', key = '--glob_outcorr--', disabled = False, enable_events = False, default = params_old['global']['outcorr'])],
                    ]
                    ,title = 'Global Parameters', size = (200,220)),
                sg.Column(layout = [
                    [sg.Button('OK', size = (10,1), key = '--accept--', button_color = 'darkgreen')],
                    [sg.Button('Default', size = (10,1), key = '--default--', button_color = 'grey')],
                    [sg.Button('Cancel', size = (10,1), key = '--cancel--', button_color = 'darkred')]]
                    ,size = (100,250))]
                ]
    
    window_PPMS = sg.Window('Parameters for Pre-Processing', layout_PPMS, size = (800,250))
    
    while True:
        event_PPMS, values_PPMS = window_PPMS.read()
        
        if event_PPMS == '--class_scale1--':
            window_PPMS['--class_scale1_how--'].Update(disabled = not values_PPMS['--class_scale1--'])
        if event_PPMS == '--vis_scale1--':
            window_PPMS['--vis_scale1_how--'].Update(disabled = not values_PPMS['--vis_scale1--'])
        
        if event_PPMS == '--class_median--':
            window_PPMS['--class_median--'].Update(value = True)
            window_PPMS['--class_concat--'].Update(value = False)
            window_PPMS['--class_separate--'].Update(value = False)
        if event_PPMS == '--class_concat--':
            window_PPMS['--class_median--'].Update(value = False)
            window_PPMS['--class_concat--'].Update(value = True)
            window_PPMS['--class_separate--'].Update(value = False)
        if event_PPMS == '--class_separate--':
            window_PPMS['--class_median--'].Update(value = False)
            window_PPMS['--class_concat--'].Update(value = False)
            window_PPMS['--class_separate--'].Update(value = True)
        if event_PPMS == '--vis_median--':
            window_PPMS['--vis_median--'].Update(value = True)
            window_PPMS['--vis_concat--'].Update(value = False)
            window_PPMS['--vis_separate--'].Update(value = False)
            window_PPMS['--glob_outcorr--'].Update(disabled = False)
        if event_PPMS == '--vis_concat--':
            window_PPMS['--vis_median--'].Update(value = False)
            window_PPMS['--vis_concat--'].Update(value = True)
            window_PPMS['--vis_separate--'].Update(value = False)
            window_PPMS['--glob_outcorr--'].Update(disabled = True, value = False)
        if event_PPMS == '--vis_separate--':
            window_PPMS['--vis_median--'].Update(value = False)
            window_PPMS['--vis_concat--'].Update(value = False)
            window_PPMS['--vis_separate--'].Update(value = True)
            window_PPMS['--glob_outcorr--'].Update(disabled = True, value = False)
        
        if event_PPMS == '--class_scale2--':
            window_PPMS['--class_scale2_how--'].Update(disabled = not values_PPMS['--class_scale2--'])
        if event_PPMS == '--vis_scale2--':
            window_PPMS['--vis_scale2_how--'].Update(disabled = not values_PPMS['--vis_scale2--'])
        
        if event_PPMS == '--glob_missing--':
            window_PPMS['--glob_minval--'].Update(disabled = not values_PPMS['--glob_missing--'])
        if event_PPMS == '--glob_minrep--':
            window_PPMS['--glob_mincount--'].Update(disabled = not values_PPMS['--glob_minrep--'])
        
        if event_PPMS == '--accept--':
            params = action.accept(values_PPMS)
            window_PPMS.close()
            break
        if event_PPMS == '--default--':
            action.reset(window_PPMS, fract_default())
            params = fract_default()
        if event_PPMS == '--cancel--' or event_PPMS == sg.WIN_CLOSED:
            params = params_old
            window_PPMS.close()
            break
    #window_PPMS.close()
    return params





def PPMS_exec_marker (params_old):
    sg.theme('DarkTeal11')
    
    if params_old['how'] == 'exclude':
        how_exclude = True
        how_majority = False
    elif params_old['how'] == 'majority':
        how_exclude = False
        how_majority = True
    
    if params_old['what'] == 'unite':
        what_unite = True
        what_intersect = False
    if params_old['what'] == 'intersect':
        what_unite = False
        what_intersect = True
    
    
    
    marker_params = copy.deepcopy(params_old)
    
    
    
    layout_PMM = [
        [sg.Text('discrepancies:\t'),
         sg.Checkbox('exclude\t', key = '--PMM_exclude--', default = how_exclude, enable_events = True),
         sg.Checkbox('take majority', key = '--PMM_majority--', default = how_majority, enable_events = True)],
        [sg.Text('selection:\t'),
         sg.Checkbox('unite\t', key = '--PMM_unite--', default = what_unite, enable_events = True),
         sg.Checkbox('intersect', key = '--PMM_intersect--', default = what_intersect, enable_events = True)],
        [sg.Button('OK', key = '--PMM_accept--', disabled = False, enable_events = True, button_color = 'dark green'),
         sg.Button('Cancel', key = '--PMM_cancel--', disabled = False, enable_events = True, button_color = 'black')]
        ]
    
    window_PMM = sg.Window('Marker Parameters', layout_PMM, size = (400,100))
    
    while True:
        event_PMM, values_PMM = window_PMM.read()
        
        if event_PMM == '--PMM_exclude--':
            marker_params['how'] = 'esclude'
            window_PMM['--PMM_exclude--'].Update(value = True)
            window_PMM['--PMM_majority--'].Update(value = False)
        if event_PMM == '--PMM_majority--':
            marker_params['how'] = 'majority'
            window_PMM['--PMM_exclude--'].Update(value = False)
            window_PMM['--PMM_majority--'].Update(value = True)
        if event_PMM == '--PMM_unite--':
            marker_params['what'] = 'unite'
            window_PMM['--PMM_unite--'].Update(value = True)
            window_PMM['--PMM_intersect--'].Update(value = False)
        if event_PMM == '--PMM_intersect--':
            marker_params['what'] = 'intersect'
            window_PMM['--PMM_unite--'].Update(value = False)
            window_PMM['--PMM_intersect--'].Update(value = True)
                
        if event_PMM == sg.WIN_CLOSED or event_PMM == '--PMM_cancel--':
            params = copy.deepcopy(params_old)
            window_PMM.close()
            break
        if event_PMM == '--PMM_accept--':
            params = copy.deepcopy(marker_params)
            window_PMM.close()
        
    #window_PMM.close()
    return params






def PPMS_exec_TP (params_old):
    sg.theme('DarkTeal11')
    tp_params = copy.deepcopy(params_old)
    
    if params_old['imputation'] == 'normal':
        is_normal = True
    else:
        is_normal = False
    if params_old['imputation'] == 'constant':
        is_constant = True
    else:
        is_constant = False
    
    layout_TPPM = [
        [sg.Text('found in at least'),
         sg.Spin(values=(list(range(1,10))), size = (10,2), key = '--tp_mincount--', disabled = False, enable_events = False, readonly = True, initial_value = params_old['minrep'], text_color = 'black')],
        [sg.Text('Imputation of Missing Values:'),
         sg.Checkbox('by normal distribution', key = '--normal--', enable_events = True, disabled = False, default = is_normal),
         sg.Checkbox('by constant (0)', key = '--constant--', enable_events = True, disabled = False, default = is_constant)],
        [sg.Button('Accept', button_color = 'dark green', key = '--TPPM_accept--', enable_events = True, disabled = False),
        sg.Button('Cancel', key = '--TPPM_cancel--', disabled = False, enable_events = True, button_color = 'black')]
        ]
    
    window_TPPM = sg.Window('TP Parameters', layout_TPPM, size = (500,100))
    
    while True:
        event_TPPM, values_TPPM = window_TPPM.read()
        
        
        if event_TPPM == '--normal--':
            tp_params['imputation'] = 'normal'
            window_TPPM['--normal--'].Update(value = True)
            window_TPPM['--constant--'].Update(value = False)
        if event_TPPM == '--constant--':
            tp_params['imputation'] = 'constant'
            window_TPPM['--constant--'].Update(value = True)
            window_TPPM['--normal--'].Update(value = False)
        
        
        if event_TPPM == sg.WIN_CLOSED or event_TPPM == '--TPPM_cancel--':
            params = copy.deepcopy(params_old)
            window_TPPM.close()
            break
        
        if event_TPPM == '--TPPM_accept--':
            params = copy.deepcopy(tp_params)
            window_TPPM.close()
            break
        
    return params






def PPMS_exec_NN (params_old):
    
    sg.theme('DarkTeal11')
    NN_params = copy.deepcopy(params_old)
    
    if 'adam' in params_old['optimizers']:
        isadam = True
    else:
        isadam = False
    if 'rmsprop' in params_old['optimizers']:
        isrms = True
    else:
        isrms = False
    if 'sgd' in params_old['optimizers']:
        issgd = True
    else:
        issgd = False
    if params_old['class_loss'] == 'mean_squared_error':
        loss_default = 'mean squared error'
    elif params_old['class_loss'] == 'binary_crossentropy':
        loss_default = 'binary cross-entropy'
    
    
    layout_NNP = [
        [sg.Checkbox('upsampling\t\t', default = params_old['upsampling'], key = '--upsampling--', disabled = False, enable_events = True),
         sg.Combo(['noised', 'average', 'noisedaverage'], default_value = params_old['upsampling_method'], key = '--upsampling_method--', size = (23,1), disabled = not params_old['upsampling'], enable_events = True, readonly = True)],
        [sg.Text('\t\t\t\tnoise:\t'),
         sg.Spin(list(range(0,6,1)), initial_value = params_old['upsampling_noise'], size = (10,1), key = '--noise--', disabled = not params_old['upsampling'], enable_events = True, readonly = True, text_color = 'black'),
         sg.Text('* std')],
        
        [sg.Checkbox('SVM-filter', default = params_old['svm_filter'], key = '--svm_filter--', disabled = False, enable_events = True)],
        
        [sg.Text('dividing-step for profile mixing:\t'),
         sg.Spin(['-',2,4,5,10], initial_value = params_old['mixed_part'], size = (15,1), key = '--mixstep--', disabled = False, enable_events = True, readonly = True, text_color = 'black')],
        
        [sg.HSep()],
        
        [sg.Text('NN optimization method:\t\t'),
         sg.Combo(['long', 'short'], default_value = params_old['NN_optimization'], size = (25,1), key = '--optimization_method--', disabled = False, enable_events = True, readonly = True)],
        
        [sg.Text('dense layer activation:\t\t'),
         sg.Combo(['relu', 'leakyrelu'], default_value = params_old['NN_activation'], size = (25,1), key = '--dense_activation--', disabled = False, enable_events = True, readonly = True)],
        
        [sg.Text('output activation:\t\t\t'),
         sg.Combo(['linear', 'sigmoid', 'softmax'], default_value = params_old['class_activation'], size = (25,1), key = '--out_activation--', disabled = False, enable_events = True, readonly = True)],
        
        [sg.Text('loss function:\t\t\t'),
         sg.Combo(['mean squared error', 'binary cross-entropy'], default_value = loss_default, size = (25,1), key = '--loss--', disabled = False, enable_events = True, readonly = True)],
        
        [sg.Text('optimizers:\t\t\t'),
         sg.Checkbox('adam', default = isadam, key = '--adam--', disabled = False, enable_events = True),
         sg.Checkbox('rmsprop', default = isrms, key = '--rmsprop--', disabled = False, enable_events = True),
         sg.Checkbox('sgd', default = issgd, key = '--sgd--', disabled = False, enable_events = True)],
        
        [sg.Text('max. training epochs:\t\t\t'),
         sg.Spin(list(range(10,110,10)), initial_value = params_old['NN_epochs'], size = (15,1), key = '--epochs--', disabled = False, enable_events = True, readonly = True, text_color = 'black')],
        
        [sg.Text('tuning runs:\t\t\t'),
         sg.Spin(list(range(1,11,1)), initial_value = params_old['rounds'], size = (15,1), key = '--rounds--', disabled = False, enable_events = True, readonly = True, text_color = 'black')],
        
        [sg.Text('ensemble runs:\t\t\t'),
         sg.Spin(list(range(5,55,5)), initial_value = params_old['subrounds'], size = (15,1), key = '--subrounds--', disabled = False, enable_events = True, readonly = True, text_color = 'black')],
        
        [sg.HSep()],
        
        [sg.Text('TP filter cut-off:\t\t\t'),
         sg.Spin(list(range(0,100,1)), initial_value = params_old['reliability'], size = (15,1), key = '--filter_threshold--', disabled = False, enable_events = True, readonly = True, text_color = 'black')],
        
        [sg.HSep()],
        
        [sg.Button('Accept', key = '--NNP_accept--', disabled = False, enable_events = True, button_color = 'dark green'),
         sg.Button('Cancel', key = '--NNP_cancel--', disabled = False, enable_events = True, button_color = 'black')]
        
        ]
    
    window_NNP = sg.Window('NN Parameters', layout_NNP, size = (470,420))
    
    while True:
        event_NNP, values_NNP = window_NNP.read()
        
        if event_NNP == '--upsampling--':
            NN_params['upsampling'] = values_NNP['--upsampling--']
            window_NNP['--upsampling_method--'].Update(disabled = not values_NNP['--upsampling--'])
            window_NNP['--noise--'].Update(disabled = not values_NNP['--upsampling--'])
        
        if event_NNP == '--upsampling_method--':
            NN_params['upsampling_method'] = values_NNP['--upsampling_method--']
        
        if event_NNP == '--noise--':
            NN_params['upsampling_noise'] = values_NNP['--noise--']
        
        if event_NNP == '--svm_filter--':
            NN_params['svm_filter'] = values_NNP['--svm_filter--']
        
        if event_NNP == '--mixstep--':
            NN_params['mixed_part'] = values_NNP['--mixstep--']
        
        if event_NNP == '--optimization_method--':
            NN_params['NN_optimization'] = values_NNP['--optimization_method--']
        
        if event_NNP == '--dense_activation--':
            NN_params['NN_activation'] = values_NNP['--dense_activation--']
        
        if event_NNP == '--out_activation--':
            NN_params['class_activation'] = values_NNP['--out_activation--']
        
        if event_NNP == '--loss--':
            if values_NNP['--loss--'] == 'mean squared error':
                NN_params['class_loss'] = 'mean_squared_error'
            elif values_NNP['--loss--'] == 'binary cross-entropy':
                NN_params['class_loss'] = 'binary_crossentropy'
        
        if event_NNP == '--adam--':
            if values_NNP['--adam--'] == True:
                NN_params['optimizers'].append('adam')
            elif values_NNP['--adam--'] == False:
                NN_params['optimizers'].remove('adam')
        if event_NNP == '--rmsprop--':
            if values_NNP['--rmsprop--'] == True:
                NN_params['optimizers'].append('rmsprop')
            elif values_NNP['--rmsprop--'] == False:
                NN_params['optimizers'].remove('rmsprop')
        if event_NNP == '--sgd--':
            if values_NNP['--sgd--'] == True:
                NN_params['optimizers'].append('sgd')
            elif values_NNP['--sgd--'] == False:
                NN_params['optimizers'].remove('sgd')
        
        if event_NNP == '--epochs--':
            NN_params['NN_epochs'] = values_NNP['--epochs--']
        
        if event_NNP == '--rounds--':
            NN_params['rounds'] = values_NNP['--rounds--']
        
        if event_NNP == '--subrounds--':
            NN_params['subrounds'] = values_NNP['--subrounds--']
        
        if event_NNP == '--filter_threshold--':
            NN_params['reliability'] = values_NNP['--filter_threshold--']
        
        
        
        if event_NNP == sg.WIN_CLOSED or event_NNP == '--NNP_cancel--':
            params = copy.deepcopy(params_old)
            window_NNP.close()
            break
        
        if event_NNP == '--NNP_accept--':
            params = copy.deepcopy(NN_params)
            window_NNP.close()
            break
    
    return params




# def NN_default ():
#     params_default = {
#         'upsampling' : True,
#         'upsampling_method' : 'noisedaverage',
#         'upsampling_noise' : 2,
#         'AE' : 'none',
#         'AE_activation' : 'leakyrelu',
#         'AE_out' : 'sigmoid',
#         'AE_epochs' : 20,
#         'svm_filter' : False,
#         'mixed_part' : 4,
#         'mixed_batch' : 0.05,
#         'NN_optimization' : 'long',
#         'NN_activation' : 'relu',
#         'class_activation' : 'linear',
#         'class_loss' : 'mean_squared_error',
#         'regularization' : 'none',
#         'optimizers' : ['adam', 'rmsprop', 'sgd'],
#         'NN_epochs' : 20,
#         'rounds' : 3,
#         'subrounds' : 10,
#         'reliability' : 95
#         }
#     return params_default



# optimizer_classes = {
#         'adam' : tf.keras.optimizers.Adam,
#         'rmsprop' : tf.keras.optimizers.RMSprop,
#         'sgd' : tf.keras.optimizers.SGD
        # }


## 'upsampling_method' : ['noised', 'average', 'noisedaverage']
## 'upsampling_noise' : [0.5, 1., 1.5, 2., 2.5]
## 'AE' : ['none', 'lite', 'full', 'full_lite']
## 'AE_activation' : ['relu', 'leakyrelu']
## 'AE_out' : ['sigmoid', relu', 'softmax', 'leakyrelu']
## 'svm_filter' : [True, False]
## 'mixed_part' : ['none', 2, 4, 5, 10]
## 'mixed_batch' : [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
## 'NN_optimization' : ['short', 'long']
## 'NN_activation' : ['relu', 'leakyrelu']
## 'class_activation' : ['sigmoid', 'softmax', 'linear']
## 'class_loss' : ['binary_crossentropy', 'mean_squared_error']
## 'regularization' : ['none', 'l1', 'l2', 'elastic']
## 'optimizers' : ['adam', 'rmsprop', 'sgd']
## 'NN_epochs' : [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]



