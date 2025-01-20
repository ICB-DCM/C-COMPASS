"""Preparameters"""


def fract_default():
    params_default = {
        "class": {
            "scale1": [
                True,
                "area",
            ],
            "corrfilter": False,
            "scale2": [False, "area"],
            "zeros": True,
            "combination": "separate",
        },
        "vis": {
            "scale1": [
                True,
                "minmax",
            ],
            "corrfilter": False,
            "scale2": [True, "minmax"],
            "zeros": True,
            "combination": "median",
        },
        "global": {
            "missing": [True, "1"],
            "minrep": [True, "2"],
            "outcorr": False,
        },
    }
    return params_default


def tp_default():
    params_default = {"minrep": 2, "imputation": "normal"}
    return params_default


def NN_default():
    params_default = {
        "upsampling": True,
        "upsampling_method": "noisedaverage",
        "upsampling_noise": 2,
        "AE": "none",
        "AE_activation": "leakyrelu",
        "AE_out": "sigmoid",
        "AE_epochs": 20,
        "svm_filter": False,
        "mixed_part": 4,
        "mixed_batch": 0.05,
        "NN_optimization": "long",
        "NN_activation": "relu",
        "class_activation": "linear",
        "class_loss": "mean_squared_error",
        "regularization": "none",
        "optimizers": ["adam", "rmsprop", "sgd"],
        "NN_epochs": 20,
        "rounds": 3,
        "subrounds": 10,
        "reliability": 95,
    }
    return params_default


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
