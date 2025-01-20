"""Core classes and functions for the ccompass package."""

from typing import Literal

from pydantic import BaseModel


class NeuralNetworkParametersModel(BaseModel):
    """Hyperparameters for the neural network."""

    upsampling: bool = True
    upsampling_method: Literal["noised", "average", "noisedaverage"] = (
        "noisedaverage"
    )
    upsampling_noise: int = 2
    AE: Literal["none", "lite", "full", "full_lite"] = "none"
    AE_activation: Literal["relu", "leakyrelu"] = "leakyrelu"
    AE_out: Literal["sigmoid", "relu", "softmax", "leakyrelu"] = "sigmoid"
    AE_epochs: int = 20
    svm_filter: bool = False
    # can be "none"
    mixed_part: int | str = 4
    mixed_batch: float = 0.05
    NN_optimization: Literal["short", "long"] = "long"
    NN_activation: Literal["relu", "leakyrelu"] = "relu"
    class_activation: Literal["sigmoid", "softmax", "linear"] = "linear"
    class_loss: Literal["binary_crossentropy", "mean_squared_error"] = (
        "mean_squared_error"
    )
    regularization: Literal["none", "l1", "l2", "elastic"] = "none"
    optimizers: list[Literal["adam", "rmsprop", "sgd"]] = [
        "adam",
        "rmsprop",
        "sgd",
    ]
    NN_epochs: int = 20
    rounds: int = 3
    subrounds: int = 10
    reliability: int = 95
