#! /usr/bin/python3

"""
@author: Vicente Yáñez
@date: Primer semestre 2016
Centro Sismologico Nacional
Universidad de Chile
Santiago, Chile.

Modificado de script de Matlab de Mike Bevis v1.0 3 Agosto 2010
Reorganizado para ser usado en servidor andes_3db
"""

import pdb

import numpy as np
from scipy.optimize import curve_fit

from gfa.gnss_analysis.hvsd import hvsd


class ModeloTrayectoria():
    """
    Clase que implementa modelo standar linear trajectory model (SLTM) y
    extended linear trajectory model (ELTM) de Bevis y Brown (2014)
    para series de tiempo.

    Para calcular los saltos requiere la funcion Heavyside(archivo hvsd.py)

    Contiene metodos para cada uno de los componentes,
        trend   : polinomio en (t-tref) de grado np
        jump    : serie Heavyside que representa saltos o discontinuidades
                  en las series de tiempo
        cycle   : serie de Fourier truncada que contiene 2*nf terminos que
                  representan oscilaciones (generalmente anuales)
        lgt     : logaritmo transitorio nl de la forma log(1 + del_t/tsc) donde
                  tsc es la escala de tiempo transitoria, y del_t=(t - t0) d
                  onde t0 es el tiempo inicial del periodo transitorio.

    La clase permite construir un array con (n+1)+nj+2*nf+nl columnas.
    El constructor trae valores por defecto para los parametros.
    Para más detalle revisar constructor.

    Construccion:
        t   : Vector de tiempo de longitud m

    Parametros [] ---> se omite componente // Input en numpy array
        tref    : tiempo de referencia. Por defecto min(t)
        n       : grado del polinomio. Debe ser un entero. Por defecto 1.
                  vacio se indica con []
        tjmp    : vector de longitud nj que contiene el tiempo de los saltos.
                  Por defecto [].
        fperiods: vector de longitud nf que contiene los periodos de la serie
                  de Fourier. Generalmente fperiods(1) = 1 año. Por defecto []
        tlt     : vector de longitud nl contiene el/los tiempo(s) de inicio de
                  la curva logaritmica. Por defecto []
        tsc     : vector de longitud nl que contiene la escala de tiempo
                  para la curva logaritmica.
                  Generalmente buenos resultados con 1 constante
                  Por defecto [1]
    """
    def __init__(self, t):
        # Parametros por defecto
        self.t = t
        self.tref = self.t.min()
        self.n = 1
        self.tjump = np.array([])
        self.fperiods = np.array([])
        self.tsc = np.array([1])
        self.tlt = np.array([])

        self.resultado = []

        # matrices vacias
        self.AT = np.array([])
        self.AF = np.array([])
        self.AJ = np.array([])
        self.AL = np.array([])

    def trend(self):
        """ Crea el array AT que es el componente trend del modelo de
        trayectoria.
        n = 0 ; modelo de posicion constante y = a
        n = 1 : modelo de velocidad constante y = a + b*dt
        n = 2 ; modelo de aceleracion constante y = a + b*dt + c*dt^2
        etc...
        donde dt = (t-tref)
        """
        # verifica que fperiods no este vacio
        if not self.n:
            print("No se ha asignado n. Modelo de trajectoria vacio")
            return self.AT
        else:
            if type(self.n) != int and self.n < 0 or self.n > 6:
                print("Error, n debe ser un entero entre 0 y 6")
                return
            m = len(self.t)
            self.AT = np.ones((m, self.n+1))
            dt = self.t - self.tref
            if self.n > 0:
                for i in range(self.n):
                    self.AT[:, i] = dt**(i+1)

    def jump(self):
        """
        Crea un array AJ con el/los componentes de "salto" del modelo de
        trayectoria.
        Necesita que tjmp no este vacio
        """
        if not self.tjump.any():
            print("No se han incorporados saltos")
            return self.AJ
        else:
            n = len(self.t)
            m = len(self.tjump)
            self.AJ = np.zeros((n, m))  # array de 0 de longitud n*m
            for i in range(m):
                htt = hvsd(self.t, self.tjump[i])
                self.AJ[:, i] = htt

    def cycle(self):
        """Crea un array AF con el/los componentes oscilatorios del modelo
        de trayectoria. Consiste en 2*nf terminos.

        Necesita que antes este establecido fperiods
        """
        # verifica que fperiods no este vacio
        if not self.fperiods.any():
            print("No se ha entregado periodo de oscilaciones")
            return self.AF
        else:
            m = len(self.t)
            # num de peridodos usadas para modelar las oscilaciones
            nf = len(self.fperiods)
            # w = array con vel. angulares para cada periodo
            w = 2*np.pi/self.fperiods
            self.AF = np.ones((m, 2*nf))  # array de 1
            # construccion AF(mx2nf)
            for i in range(nf):
                self.AF[:, i-1] = np.sin(w[i]*self.t)
                self.AF[:, i] = np.cos(w[i]*self.t)

    def lgt(self):
        """
        Crea un array AL que contiene los componentes para la parte logaritmica
        del modelo de trayectoria, asumiendo que la escala de tiempo de la
        curva logaritmica es ya conocida.

        Necesita las variables de clase tlt y tsc
        """
        if not self.tlt.any() or not self.tsc.any():
            print("tlt o tsc esta vacio")
            return
        else:
            nl = len(self.tlt)
            self.AL = np.zeros((len(self.t), nl))
            for i in range(nl):
                ib = np.where(self.t <= self.tlt[i])
                ia = np.where(self.t > self.tlt[i])
                self.AL[ib, i] = np.zeros((len(ib), 1))
                self.AL[ia, i] = np.log10(1 + (self.t[ia]-self.tlt[i]) /
                                          self.tsc[i])

    def modelo_trayectoria(self, y):
        """
        Metodo que se encarga de solucionar el modelo de trayectorias
        usando solucion de los minimos cuadrados.
        Primero construye el array len(t) x num_parametros que es necesario
        para llamar a la funcion que resuelve los minimos cuadrados.

        """
        modelo = []
        # verifica que no este vacio
        if not self.n:
            print("No se ha asignado n. Modelo de trayectoria vacio")
        else:
            self.trend()
            modelo = self.AT
        if not self.tjump.any():
            print("No se han incorporados saltos")
        else:
            self.jump()
            # en caso de que no se haya asignado un valor previamente
            if not modelo.any():
                modelo = self.AJ
            else:
                modelo = np.hstack((modelo, self.AJ))
        if not self.fperiods.any():
            print("No se ha entregado periodo de oscilaciones")
        else:
            self.cycle()
            if not modelo.any():
                modelo = self.AF
            else:
                modelo = np.hstack((modelo, self.AF))
        if not self.tlt.any() or not self.tsc.any():
            print("tlt o tsc esta vacio")
        else:
            self.lgt()
            if not modelo.any():
                modelo = self.AL
            else:
                modelo = np.hstack((modelo, self.AL))

        # numero de parametros por componente
        nt = self.n + 1
        nj = len(self.tjump)
        nf = 2*len(self.fperiods)
        nl = len(self.tlt)

        # posiciones en array
        n2 = nt+nj
        n3 = nt+nj+nf
        n4 = nt+nj+nf+nl

        # solucion mediante minimos cuadrados
        parametros, residual, rank, s = np.linalg.lstsq(modelo, y)

        # devuelve el valor de distancia del modelo para cada componente
        res_nt = np.matrix(modelo[:, :nt])*np.matrix(parametros[:nt]).T
        res_nj = np.matrix(modelo[:, nt:n2])*np.matrix(parametros[nt:n2]).T
        res_nf = np.matrix(modelo[:, n2:n3])*np.matrix(parametros[n2:n3]).T
        res_nl = np.matrix(modelo[:, n3:n4])*np.matrix(parametros[n3:n4]).T
        # suma cada componente para el resultado final
        res_nt = np.array(res_nt)
        res_nj = np.array(res_nj)
        res_nf = np.array(res_nf)
        res_nl = np.array(res_nl)

        self.resultado = [res_nt, res_nj, res_nf, res_nl]
        self.total = res_nt.T[0] + res_nj.T[0] + res_nf.T[0] + res_nl.T[0]

        return self.resultado, residual

    def save_components(self, dir_mod, estac, eje):
        """
        Guarda el modelo obtenido, requiere haber construido el modelo
        previamente

        Input
        dir_mod     :   path del Modelo str
        estac       :   nombre estacion str
        eje         :   eje de la estacion que se va a guardar str

        Guarda archivos con dos columnas: tiempo  desplazamiento
        """
        array_guardado = []
        nt_guardar = []
        nf_guardar = []
        nl_guardar = []
        nj_guardar = []

        # formatear array para guardado
        for i in range(len(self.t)):
            nt_guardar.append([self.t[i], self.resultado[0][i]])
            nj_guardar.append([self.t[i], self.resultado[1][i]])
            nf_guardar.append([self.t[i], self.resultado[2][i]])
            nl_guardar.append([self.t[i], self.resultado[3][i]])

        # desplazamiento modelo por partes
        np.savetxt(dir_mod+'comp/'+estac+eje+'_nt.txt', nt_guardar, '%12.8f')
        np.savetxt(dir_mod+'comp/'+estac+eje+'_nj.txt', nj_guardar, '%12.8f')
        np.savetxt(dir_mod+'comp/'+estac+eje+'_nf.txt', nf_guardar, '%12.8f')
        np.savetxt(dir_mod+'comp/'+estac+eje+'_nl.txt', nl_guardar, '%12.8f')
        return
