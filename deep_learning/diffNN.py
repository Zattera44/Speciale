import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt

class ForwardLayer(tf.keras.layers.Layer):
    def __init__(self, units, **kwargs):
        super().__init__(**kwargs)
        self.units = units
        
    def build(self, input_shape):
        self.w = self.add_weight(
            name = 'weights',
            shape=(input_shape[-1], self.units),
            initializer="glorot_normal",
            trainable=True
        )
        self.b = self.add_weight(
            name = 'bias',
            shape=(self.units,), 
            initializer="zeros", 
            trainable=True
        )
        super().build(input_shape)

    def call(self, inputs, first=False):
        if not first:
            return tf.keras.activations.softplus(inputs) @ self.w + self.b
        if first:
            return inputs @ self.w + self.b

class BackpropLayer(tf.keras.layers.Layer):
    def __init__(self, twin: ForwardLayer, units=20, **kwargs):
        super().__init__(**kwargs)
        self.twin = twin
        self.units = self.twin.units
        self.units = units

    def build(self, input_shape):
        self.w = self.twin.w
        super().build(input_shape)

    def call(self, inputs, output=False, first=False):
        if first:
            return tf.transpose(self.w) * tf.keras.activations.sigmoid(inputs)
        if not output and not first:
            z, zbar = inputs
            return zbar @ tf.transpose(self.w) * tf.keras.activations.sigmoid(z)
        if output and not first:
            return inputs @ tf.transpose(self.w)

def create_graph(layers=4, units=20):
    forward_layers, backprop_layers, z, zbar = {}, {}, {}, {}
    for l in range(layers):
        name = 'forward_{}'.format(l+1)
        new_layer = ForwardLayer(units=units, name=name)
        forward_layers[name] = new_layer
    forward_layers['output'] = ForwardLayer(units=1, name='output')
    l = 0
    for key in reversed(forward_layers):
        if key == 'forward_1':
            break
        else: 
            name = 'backprop_{}'.format(l+1)
            new_layer = BackpropLayer(twin = forward_layers[key], name=name)
            backprop_layers[name] = new_layer
            l += 1
    backprop_layers['derivs'] = BackpropLayer(twin = forward_layers[key], name='derivs')

    z['z_0'] = tf.keras.layers.Input(shape=(1,), name='input')
    old_name = 'z_0'
    l = 1
    for key in forward_layers:
        layer = forward_layers[key]
        if key == 'forward_1':
            first = True
        else:
            first = False
        if key =='output':
            y = layer(z[new_name])
        else:
            new_name = 'z_{}'.format(l)
            z[new_name] = layer(z[old_name], first=first)
            old_name = new_name
            l += 1
    l = layers
    for layer_key, z_key in zip(backprop_layers, reversed(z)):
        layer = backprop_layers[layer_key]
        if l == 0:
            derivs = layer(inputs = zbar[old_name], output=True)
        new_name = 'zbar_{}'.format(l)
        if layer_key == 'backprop_1':
            zbar[new_name] = layer(z[z_key], first = True)
            old_name = new_name
        else:
            zbar[new_name] = layer([z[z_key],zbar[old_name]])
            old_name = new_name
        l -= 1
    return tf.keras.Model(inputs=z['z_0'], outputs=[y,derivs])  

def normalize_data(x_raw, y_raw, dydx_raw, epsilon = 1.0e-08):
    x_mean = x_raw.mean(axis=0)
    x_std = x_raw.std(axis=0) + epsilon
    x = (x_raw- x_mean) / x_std
    y_mean = y_raw.mean(axis=0)
    y_std = y_raw.std(axis=0) + epsilon
    y = (y_raw-y_mean) / y_std
    dydx = dydx_raw / y_std * x_std 
    lambda_j = 1.0 / np.sqrt((dydx ** 2).mean(axis=0)).reshape(1, -1)
    return x_mean, x_std, x, y_mean, y_std, y, dydx, lambda_j 

def create_derivsLoss(lambda_j):
    def derivsLoss(dydx_true, dydx_pred):
        dydx_true = dydx_true  * lambda_j
        dydx_pred = dydx_pred * lambda_j
        error = dydx_true - dydx_pred
        return tf.reduce_mean(tf.square(error))
    return derivsLoss

class LRSchedule(tf.keras.optimizers.schedules.LearningRateSchedule):

  def __init__(self, n_epochs, size, batch_size):
    self.learning_rate_schedule = [   (0.0, 1.0e-8), \
                                      (0.2, 0.1),    \
                                      (0.6, 0.01),   \
                                      (0.9, 1.0e-6), \
                                      (1.0, 1.0e-8)  ]
    self.n_epochs = n_epochs
    self.lr_schedule_epochs, self.lr_schedule_rates = zip(*self.learning_rate_schedule)
    self.batch_size = batch_size
    self.size = size
    self.n_steps = self.batch_size * self.size
    self.epoch = 1
    self.counter = 0
    self.number = self.size/self.batch_size

  def __call__(self, step):
    self.counter += 1
    if self.counter % self.number == 0:
        self.epoch += 1
    lr = np.interp(self.epoch / self.n_epochs, self.lr_schedule_epochs, self.lr_schedule_rates)

    return lr  

default_schedule = [(0.0, 1.0e-8), \
                    (0.2, 0.1),    \
                    (0.6, 0.01),   \
                    (0.9, 1.0e-6), \
                    (1.0, 1.0e-8)  ]

def lr_scheduler(epoch, n_epochs, lr_schedule = default_schedule):
    lr_schedule_epochs, lr_schedule_rates = zip(*lr_schedule)
    return np.interp(epoch / n_epochs, lr_schedule_epochs, lr_schedule_rates)

def differential_mse(y_batch, yhat_batch, dydx_batch, dydxhat_batch, lambd):
    error = y_batch - yhat_batch
    error_bar = dydx_batch - dydxhat_batch
    return tf.reduce_mean(tf.square(error)) + lambd * tf.reduce_mean(tf.square(error_bar))

# pinligt dårligt skrevet funktion, men det gør hvad den skal
def plot_value_delta(xTest, yPred, yTest, dydxPred, dydxTest, size):
    fig, ax = plt.subplots(1, 2, squeeze=False, dpi=90)
    fig.set_size_inches(9.5, 4)
    ax[0,0].plot(xTest,dydxPred, 'co', markersize=2, color='red', label='Predicted')
    ax[0,0].plot(xTest,dydxTest, color='blue', label='Monte Carlo')
    ax[0,1].plot(xTest,yPred, 'co', markersize=2, color='red', label='Predicted')
    ax[0,1].plot(xTest,yTest, color='blue', label='Monte Carlo')
    ax[0,0].set_ylabel("Delta")
    ax[0,1].set_ylabel("Price")
    ax[0,0].set_xlabel("Spot")
    ax[0,1].set_xlabel("Spot")
    errors = 100*(dydxPred - dydxTest)
    rmse = np.sqrt((errors ** 2).mean(axis=0))
    t = "RMSE = %.2f" % rmse
    ax[0,0].set_title(t)
    errors = (yPred - yTest)
    rmse = np.sqrt((errors ** 2).mean(axis=0)); rmse
    t = "RMSE = %.2f" % rmse
    ax[0,1].set_title(t)
    t = "Size = %.0f" % size
    plt.suptitle(t)
    ax[0,0].legend()
    plt.show()
