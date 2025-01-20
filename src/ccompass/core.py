"""Core classes and functions for the ccompass package."""

from typing import Literal

from pydantic import BaseModel


class NeuralNetworkParametersModel(BaseModel):
    """Hyperparameters for the neural network."""

    #: Perform upsampling?
    upsampling: bool = True
    #: Method for upsampling
    upsampling_method: Literal[
        "none", "noised", "average", "noisedaverage"
    ] = "noisedaverage"
    #: Noise level for upsampling (standard deviations)
    upsampling_noise: float = 2
    # FIXME: unused
    AE: Literal["none", "lite", "full", "full_lite"] = "none"
    # FIXME: unused
    AE_activation: Literal["relu", "leakyrelu"] = "leakyrelu"
    # FIXME: unused
    AE_out: Literal["sigmoid", "relu", "softmax", "leakyrelu"] = "sigmoid"
    # FIXME: unused
    AE_epochs: int = 20
    #: Perform SVM filtering?
    svm_filter: bool = False
    #: ...
    # FIXME: can be "none"; == 0?!
    mixed_part: int | str = 4
    #: ...
    mixed_batch: float = 0.05
    #: Long or short optimization?
    NN_optimization: Literal["short", "long"] = "long"
    #: Neural network activation function
    NN_activation: Literal["relu", "leakyrelu"] = "relu"
    #: Neural network class layer activation function
    class_activation: Literal["sigmoid", "softmax", "linear"] = "linear"
    #: Neural network training loss function
    class_loss: Literal["binary_crossentropy", "mean_squared_error"] = (
        "mean_squared_error"
    )
    #: FIXME: unused
    regularization: Literal["none", "l1", "l2", "elastic"] = "none"
    #: Optimizers to include in the hyperparameter search
    optimizers: list[Literal["adam", "rmsprop", "sgd"]] = [
        "adam",
        "rmsprop",
        "sgd",
    ]
    #: Number of epochs for the neural network training
    NN_epochs: int = 20
    #: ...
    rounds: int = 3
    #: ...
    subrounds: int = 10
    #: Percentile threshold for ... ?
    reliability: int = 95
