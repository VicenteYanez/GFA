#! /usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np


def hvsd(t, t0, *tr):
    """
    Modificado y traduccion de script de matlab de Mike Bevis.
    @editor: Vicente Yáñez
    @date: Enero 2016
    Centro Sismologico Nacional
    Santiago, Chile


    Funcion heavyside.
    Hay 3 formas de uso. Implementada solo la tercera en el script.
    1. HVSD(t): uso estandar
        hvsd(t) = 0     if t < 0
        hvsd(t) = 0.5   if t = 0
        hvsd(t) = 1     if t > 0

    2. HVSD(t, t0): Evalua H(t-t0)
        hvsd(t, t0) = 0     if t < 0
        hvsd(t, t0) = 0.5   if t = 0
        hvsd(t, t0) = 1     if t > 0

    3. HVSD(t, t0, tr): evalua la expresion (H(t-t0) - H(t0-tr))
        if t0 < tr
            hvsd(t, t0, tr) = -1    if t < t0
            hvsd(t, t0, tr) = -0.5  if t = t0
            hvsd(t, t0, tr) =  0    if t > t0
        if t0 = tr
            hvsd(t, t0, tr) = -0.5   if t < t0
            hvsd(t, t0, tr) =  0     if t = t0
            hvsd(t, t0, tr) = +0.5   if t > t0
        if t0 > tr
            hvsd(t, t0, tr) = 0    if t < t0
            hvsd(t, t0, tr) = 0.5  if t = t0
            hvsd(t, t0, tr) = 1    if t > t0

    'Esta última forma es util para modelar los saltos en las series de
    tiempo (En realidad t0 es nunca igual a tr en este contexto dado que
    podemos elegir el tiempo de referencia tr, y no se elije en el
    momento en que la discontinuidad aparece)' Traducido de Bevis

    Nota: pero funcion jmp, solo usa entrega dos argumentos, por
    lo que se usa el segundo modelo.

    t0 y tr : escalares
    t       : puede ser un vector
    """

    ht = np.array(t)

    i1 = np.where(t < t0)  # posicion de valores en t que son < a t0
    if i1[0].any():
        ht[i1] = 0

    i2 = np.where(t == 0)
    if i2[0].any():
        ht[i2] = 0.5

    i3 = np.where(t > t0)
    if i3[0].any():
        ht[i3] = 1

    # corre si se entrega un valor escalar tr
    if tr:
        if t0 == tr:
            H2 = 0.5
        elif t0 > tr:
            H2 = 1.0
        else:
            H2 = 0

        ht -= H2

    return ht
