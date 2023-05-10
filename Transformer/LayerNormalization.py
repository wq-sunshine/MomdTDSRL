from keras import backend as K
from keras.engine.topology import Layer
import math
import numpy as np
import tensorflow as tf
from keras.initializers import Ones,Zeros

class LayerNormalization(Layer):
    '''
    Layer normalization of the input with trainable gain_parameter and bias parameter
    Source: https://github.com/Lsdefine/attention-is-all-you-need-keras/blob/master/transformer.py#L14
    '''
    def __init__(self, **kwargs):
        super(LayerNormalization, self).__init__(**kwargs)

    def build(self, input_shape):
        self.epsilon   = 1e-6 #大于0的小浮点数，用于防止除0错误
        self.gain_para = self.add_weight(name='gain_parameter',
                                          shape=[input_shape[1],1],
                                          initializer = Ones(),trainable=True)
        self.bias_para = self.add_weight(name='bias_parameter',
                                          shape=[input_shape[1],1],
                                          initializer = Zeros(),
                                          trainable = True)

        super(LayerNormalization, self).build(input_shape)  # Be sure to call this at the end

    def call(self, x):

        mean = K.mean(x,axis=-1,keepdims=True)
        std  = K.std(x,axis=-1,keepdims=True)#元素X的标准差

        normalised_vec = ((x-mean) / (std + self.epsilon)) * self.gain_para + self.bias_para
        return normalised_vec# 标准化向量

    def compute_output_shape(self, input_shape):
        # assert isinstance(input_shape, list)
        return input_shape
