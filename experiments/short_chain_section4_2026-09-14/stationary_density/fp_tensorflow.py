#!/usr/bin/env python3
"""TensorFlow implementation of the audited reduced adjoint residual."""

from __future__ import annotations

import tensorflow as tf


def drift_tf(state: tf.Tensor, t1: float, t3: float, gamma: float = 0.1) -> tf.Tensor:
    state = tf.cast(state, tf.float64)
    i1, i2, i3, th1, th3 = tf.unstack(state, axis=-1)
    s1, c1 = tf.sin(th1), tf.cos(th1)
    s3, c3 = tf.sin(th3), tf.cos(th3)
    g = tf.cast(gamma, tf.float64)
    t1 = tf.cast(t1, tf.float64); t3 = tf.cast(t3, tf.float64)
    b1 = 4*i1*i2*s1-4*g*i1*i2*(1+c1)+4*g*t1-2*g*(i1*i1+2*i1*i3)
    b2 = -4*i2*(i1*s1+i3*s3)
    b3 = 4*i3*i2*s3-4*g*i3*i2*(1+c3)+4*g*t3-2*g*(i3*i3+2*i1*i3)
    bt1 = 2*(i2-i1)+4*(i2-i1)*c1+4*g*i2*s1-4*i3*c3
    bt3 = 2*(i2-i3)+4*(i2-i3)*c3+4*g*i2*s3-4*i1*c1
    return tf.stack((b1,b2,b3,bt1,bt3), axis=-1)


def divergence_tf(state: tf.Tensor, gamma: float = 0.1) -> tf.Tensor:
    state = tf.cast(state, tf.float64)
    i1, i2, i3, th1, th3 = tf.unstack(state, axis=-1)
    s1, c1 = tf.sin(th1), tf.cos(th1)
    s3, c3 = tf.sin(th3), tf.cos(th3)
    g = tf.cast(gamma, tf.float64)
    return (
        4*i2*s1-4*g*i2*(1+c1)-4*g*(i1+i3)
        -4*(i1*s1+i3*s3)
        +4*i2*s3-4*g*i2*(1+c3)-4*g*(i3+i1)
        -4*(i2-i1)*s1+4*g*i2*c1
        -4*(i2-i3)*s3+4*g*i2*c3
    )


def relative_adjoint_residual(model, state: tf.Tensor, t1: float, t3: float,
                              gamma: float = 0.1) -> tf.Tensor:
    state = tf.cast(state, tf.float64)
    with tf.GradientTape(persistent=True) as outer:
        outer.watch(state)
        with tf.GradientTape() as inner:
            inner.watch(state)
            log_prob = model.log_prob(state)
            log_prob_sum = tf.reduce_sum(log_prob)
        grad = inner.gradient(log_prob_sum, state)
        gradient_sums = [tf.reduce_sum(grad[:, index]) for index in range(5)]
    hdiag = []
    for index in range(5):
        derivative = outer.gradient(gradient_sums[index], state)
        hdiag.append(derivative[:, index])
    del outer
    hdiag = tf.stack(hdiag, axis=-1)
    b = drift_tf(state, t1, t3, gamma)
    divb = divergence_tf(state, gamma)
    i1, _, i3, _, _ = tf.unstack(state, axis=-1)
    g = tf.cast(gamma, tf.float64)
    t1t = tf.cast(t1, tf.float64); t3t = tf.cast(t3, tf.float64)
    a = tf.stack((
        4*g*t1t*i1, tf.zeros_like(i1), 4*g*t3t*i3,
        4*g*t1t/i1, 4*g*t3t/i3,
    ), axis=-1)
    residual = -divb-tf.reduce_sum(b*grad, axis=-1)
    residual += tf.reduce_sum(a*(hdiag+grad*grad), axis=-1)
    residual += 8*g*t1t*grad[:, 0]+8*g*t3t*grad[:, 2]
    return residual
