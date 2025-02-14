"""The neural network model for multi-compartment classification."""

import keras_tuner as kt
import tensorflow as tf

# noinspection PyUnresolvedReferences
from tensorflow import keras

# noinspection PyUnresolvedReferences
from tensorflow.keras import ops

# noinspection PyUnresolvedReferences
from tensorflow.keras.backend import epsilon

from .core import NeuralNetworkParametersModel

#: mapping of optimizer names to optimizer classes
# noinspection PyUnresolvedReferences
optimizer_classes = {
    # noinspection PyUnresolvedReferences
    "adam": tf.keras.optimizers.Adam,
    # noinspection PyUnresolvedReferences
    "rmsprop": tf.keras.optimizers.RMSprop,
    # noinspection PyUnresolvedReferences
    "sgd": tf.keras.optimizers.SGD,
}


class FNN_Classifier(kt.HyperModel):
    """The neural network hypermodel for multi-compartment classification."""

    def __init__(
        self,
        nn_params: NeuralNetworkParametersModel,
        fixed_hp=None,
        set_shapes=None,
    ):
        super().__init__()
        self.fixed_hp = fixed_hp
        self.set_shapes = set_shapes
        self.chosen_hp = {}
        self.nn_params = nn_params

    def build(self, hp):
        model = keras.Sequential()
        # Input layer, size is the number of fractions
        # noinspection PyUnresolvedReferences
        model.add(
            tf.keras.Input(
                (self.set_shapes[0],),
            )
        )

        # fixed or tunable hyperparameters
        if self.fixed_hp:
            optimizer_choice = self.fixed_hp["optimizer"]
            learning_rate = self.fixed_hp["learning_rate"]
            units = self.fixed_hp["units"]
        else:
            optimizer_choice = hp.Choice(
                "optimizer", self.nn_params.optimizers
            )
            learning_rate = hp.Float(
                "learning_rate",
                min_value=1e-4,
                max_value=1e-1,
                sampling="log",
            )
            if self.nn_params.NN_optimization == "short":
                units = hp.Int(
                    "units",
                    min_value=int(
                        min(self.set_shapes)
                        + 0.4 * (max(self.set_shapes) - min(self.set_shapes))
                    ),
                    max_value=int(
                        min(self.set_shapes)
                        + 0.6 * (max(self.set_shapes) - min(self.set_shapes))
                    ),
                    step=2,
                )
            elif self.nn_params.NN_optimization == "long":
                units = hp.Int(
                    "units",
                    min_value=min(self.set_shapes),
                    max_value=max(self.set_shapes),
                    step=2,
                )
            else:
                raise ValueError(
                    f"Unknown optimization: {self.nn_params.NN_optimization}"
                )

        # dense layer 1 with tunable size
        if self.nn_params.NN_activation == "relu":
            model.add(keras.layers.Dense(units, activation="relu"))
        elif self.nn_params.NN_activation == "leakyrelu":
            hp_alpha = hp.Float(
                "alpha", min_value=0.05, max_value=0.3, step=0.05
            )
            model.add(keras.layers.Dense(units))
            model.add(keras.layers.LeakyReLU(hp_alpha))

        # dense layer 2 with size according to the number of compartments
        model.add(
            keras.layers.Dense(
                self.set_shapes[1],
                activation=self.nn_params.class_activation,
            )
        )
        model.add(keras.layers.ReLU())

        # normalization layer
        model.add(keras.layers.Lambda(sum1_normalization))

        optimizer = optimizer_classes[optimizer_choice](
            learning_rate=learning_rate
        )
        # noinspection PyUnresolvedReferences
        model.compile(
            loss=self.nn_params.class_loss,
            optimizer=optimizer,
            metrics=[
                tf.keras.metrics.MeanSquaredError(),
                tf.keras.metrics.MeanAbsoluteError(),
            ],
        )

        if not self.fixed_hp:
            self.chosen_hp = {
                "optimizer": optimizer_choice,
                "learning_rate": learning_rate,
                "units": units,
            }

        return model

    def get_chosen_hyperparameters(self):
        return self.chosen_hp


def sum1_normalization(x):
    """Normalize the input to sum to 1."""
    return x / (ops.sum(x, axis=1, keepdims=True) + epsilon())
