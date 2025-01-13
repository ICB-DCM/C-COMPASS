import PySimpleGUI as sg




## parameters:
#------------------------------------------------------------------------------





## internal functions:
#------------------------------------------------------------------------------





## external functions:
#------------------------------------------------------------------------------

def set_combinations (params):
    combination_both = {}
    if params['class']['combination'] == 'median':
        combination_both['class_median'] = True
        combination_both['class_concat'] = False
        combination_both['class_separate'] = False
    elif params['class']['combination'] == 'concat':
        combination_both['class_median'] = False
        combination_both['class_concat'] = True
        combination_both['class_separate'] = False
    elif params['class']['combination'] == 'separate':
        combination_both['class_median'] = False
        combination_both['class_concat'] = False
        combination_both['class_separate'] = True
    if params['vis']['combination'] == 'median':
        combination_both['vis_median'] = True
        combination_both['vis_concat'] = False
        combination_both['vis_separate'] = False
    elif params['vis']['combination'] == 'concat':
        combination_both['vis_median'] = False
        combination_both['vis_concat'] = True
        combination_both['vis_separate'] = False
    elif params['vis']['combination'] == 'separate':
        combination_both['vis_median'] = False
        combination_both['vis_concat'] = False
        combination_both['vis_separate'] = True
    return (combination_both)

def accept (values):
    params = {
        'class' : {
            'scale1' : [
                values['--class_scale1--'],
                values['--class_scale1_how--'],
                ],
            'corrfilter' : values['--class_corrfilter--'],
            'scale2' : [
                values['--class_scale2--'],
                values['--class_scale2_how--']
                ],
            'zeros' : values['--class_zeros--']
            },
        'vis' : {
            'scale1' : [
                values['--vis_scale1--'],
                values['--vis_scale1_how--'],
                ],
            'corrfilter' : values['--vis_corrfilter--'],
            'scale2' : [
                values['--vis_scale2--'],
                values['--vis_scale2_how--']
                ],
            'zeros' : values['--vis_zeros--']
            },
        'global' : {
            'missing' : [
                values['--glob_missing--'],
                values['--glob_minval--']
                ],
            'minrep' : [
                values['--glob_minrep--'],
                values['--glob_mincount--']
                ],
            'outcorr' : values['--glob_outcorr--']
            }
        }
    if values['--class_median--']:
        params['class']['combination'] = 'median'
    elif values['--class_concat--']:
        params['class']['combination'] = 'concat'
    elif values['--class_separate--']:
        params['class']['combination'] = 'separate'
    if values['--vis_median--']:
        params['vis']['combination'] = 'median'
    elif values['--vis_concat--']:
        params['vis']['combination'] = 'concat'
    elif values['--vis_separate--']:
        params['vis']['combination'] = 'separate'
    return params

def reset (window, params):
    window['--class_scale1--'].Update(value = params['class']['scale1'][0])
    window['--class_scale1_how--'].Update(value = params['class']['scale1'][1])
    window['--class_corrfilter--'].Update(value = params['class']['corrfilter'])
    window['--class_scale2--'].Update(value = params['class']['scale2'][0])
    window['--class_scale2_how--'].Update(value = params['class']['scale2'][1])
    window['--class_zeros--'].Update(value = params['class']['zeros'])
    
    window['--vis_scale1--'].Update(value = params['vis']['scale1'][0])
    window['--vis_scale1_how--'].Update(value = params['vis']['scale1'][1])
    window['--vis_corrfilter--'].Update(value = params['vis']['corrfilter'])
    window['--vis_scale2--'].Update(value = params['vis']['scale2'][0])
    window['--vis_scale2_how--'].Update(value = params['vis']['scale2'][1])
    window['--vis_zeros--'].Update(value = params['vis']['zeros'])
    
    window['--glob_missing--'].Update(value = params['global']['missing'][0])
    window['--glob_minval--'].Update(value = params['global']['missing'][1])
    window['--glob_minrep--'].Update(value = params['global']['minrep'][0])
    window['--glob_mincount--'].Update(value = params['global']['minrep'][1])
    window['--glob_outcorr--'].Update(value = params['global']['outcorr'])
    
    if params['class']['combination'] == 'median':
        window['--class_median--'].Update(value = True)
        window['--class_concat--'].Update(value = False)
        window['--class_separate--'].Update(value = False)
    elif params['class']['combination'] == 'concat':
        window['--class_median--'].Update(value = False)
        window['--class_concat--'].Update(value = True)
        window['--class_separate--'].Update(value = False)
    elif params['class']['combination'] == 'separate':
        window['--class_median--'].Update(value = False)
        window['--class_concat--'].Update(value = False)
        window['--class_separate--'].Update(value = True)
    if params['vis']['combination'] == 'median':
        window['--vis_median--'].Update(value = True)
        window['--vis_concat--'].Update(value = False)
        window['--vis_separate--'].Update(value = False)
    elif params['vis']['combination'] == 'concat':
        window['--vis_median--'].Update(value = False)
        window['--vis_concat--'].Update(value = True)
        window['--vis_separate--'].Update(value = False)
    elif params['vis']['combination'] == 'separate':
        window['--vis_median--'].Update(value = False)
        window['--vis_concat--'].Update(value = False)
        window['--vis_separate--'].Update(value = True)
    return










































