#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math as m
import pandas as pd
import numpy as np
import scipy.optimize as spo
import Metrics as Mt
import matplotlib.pyplot as plt

# =============================== Definicion para cada modelo de la función parte derecha ============================ #


def fpd01(ci, x, k, tau):
    """
    :param ci: Valor de la condicion inicial de la variable de estado (Caudal)
    :param x: Entrada del sistema (Precipitacion)
    :param k: Parametro de eficiencia (Coeficiente de escorrentia)
    :param tau: Tiempo de relajamiento (Tiempo de concentración)
    :return: Caudal de salida en el tiempo
    """
    return -ci / (k * tau) + x / tau


def fpd02(ci, x, omega, a, hs, s, n):
    """"
    :param ci: Valor de la condicion inicial de la variable de estado
    :param x: Entrada del sistema
    :param omega: Ancho medio de la cuenca
    :param a: Area de la cuenca
    :param hs: Capacidad de almacenamiento
    :param s: Pendiente media de la cuenca
    :param n: Rugosidad
    :return: Caudal de salida de la cuenca
    """
    return omega / (a * n) * ((ci - hs + x - omega / (a * n) * ((ci - hs)**(5/3)) * s**(1/2))**(5/3)) * s**(1/2)


def fpd03(ci, x, k1):
    """"
    :param ci: Valor de la condicion inicial de la variable de estado
    :param x: Entrada del sistema
    :param k1: Ancho medio de la cuenca
    :return: Caudal de salida de la cuenca
    """
    return ci / (k1 * x) - (1 - ci / (k1 * x)) * m.log(1 - ci / (k1 * x))


# ================================================ Definicion de modelos ============================================= #
def euler(dt, ci, fpd):
    """
    :param dt: Discretizacion temporal
    :param ci: Condicion inicial
    :param fpd: Funcion parte derecha
    :return: Q (t+dt)
    """
    return fpd * dt + ci


def run_euler(datos, k, tau):
    dt = 1.0
    nt = datos['Q'].size
    Q = datos['Q'].values
    x = datos['X'].values
    Qeuler = np.zeros(nt)
    ci = Q[0]
    Qeuler[0] = ci

    # Ciclo de calculo en el intervalo de tiempo
    for i in range(1, nt):
        fpd = fpd01(ci, x[i], k, tau)
        Qeuler[i] = euler(dt, ci, fpd)
        ci = Q[i]

    return Qeuler


def run_euler2(datos, hs, s, n):
    dt = 1.0
    # Area de la cuenca en m2
    area = 1000000000
    omega = 20000
    nt = datos['Q'].size
    q = datos['Q'].values
    x = datos['X'].values
    Qeuler2 = np.zeros(nt)
    ci = x[0] / 1000 - (q[0] * (30 * 86400)) / area
    Qeuler2[0] = q[0]

    # Ciclo de calculo en el intervalo de tiempo
    for i in range(1, nt):
        fpd = fpd02(ci, x[i], omega, area, hs, s, n)
        Qeuler2[i] = euler(dt, ci, fpd)
        ci = x[i] / 1000 - (q[i] * (30 * 86400)) / area

    return Qeuler2

# ========================================== Ejemplos uso de modelos numericos ======================================= #
if __name__ == "__main__":
    input_file = 'Datos_Betania.xlsx'
    datos = pd.read_excel(input_file, index_col=0)

    # ************************************ Ejemplo modelo de Euler con fpd01 ************************************ #
    print 223 * '*'
    print 'Ejemplo modelo de Euler con fpd01'
    k = 0.4
    tau = 0.8
    print 'parametros iniciales'
    print 'k=', k, 'tau=', tau

    # Optimización de parametros (2)
    par_ini = [k, tau]
    Qobs = datos['Q'].values
    par_opt, par_cov = spo.curve_fit(run_euler, datos, Qobs, par_ini)
    datos['Euler_fpd01'] = run_euler(datos, *par_opt)
    print 'parametros optimizados'
    print 'k=', par_opt[0], 'tau=', par_opt[1]

    print datos.to_string()

    # Guardado de los resultados como un archivo de texto con numpy
    np.savetxt('Simu_Euler_fpd01.txt', datos.values, fmt='%.18e', delimiter=' ',
               newline='\n', header='', footer='', comments='# ')

    # Calculo de las metricas de desempeño
    sim_metricas = np.loadtxt('Simu_Euler_fpd01.txt')
    Obs = sim_metricas[:, 1:2]
    Sim = sim_metricas[:, 2:3]
    metricas_euler1 = Mt.run_metricas(Obs, Sim)
    # print 223 * '*'
    # print '///Metricas Euler con fpd01///'
    # print metricas_euler1.to_string()
    # print 223 * '*'

    # ************************************ Ejemplo modelo de Euler con fpd02 ************************************ #
    print 25 * '*', 'Ejemplo modelo de Euler con fpd02', 25 * '*'
    hs = 0.25
    s = 0.001
    n = 0.04
    print 'parametros iniciales'
    print 'hs=', hs, 's=', s, 'n=', n

    # Optimización de parametros (2)
    par_ini2 = [hs, s, n]
    Qobs = datos['Q'].values
    par_opt2, par_cov2 = spo.curve_fit(run_euler2, datos, Qobs, par_ini2)
    datos['Euler_fpd02'] = run_euler2(datos, *par_opt2)
    print 'parametros optimizados'
    print 'hs=',par_opt2[0], 's=', par_opt2[1], 'n=', par_opt2[2]

    print datos.to_string()

    # Guardado de los resultados como un archivo de texto con numpy
    np.savetxt('Simu_Euler_fpd02.txt', datos.values, fmt='%.18e', delimiter=' ',
               newline='\n', header='', footer='', comments='# ')

    # Calculo de las metricas de desempeño
    sim_metricas2 = np.loadtxt('Simu_Euler_fpd02.txt')
    Obs2 = sim_metricas2[:, 1:2]
    Sim2 = sim_metricas2[:, 3:4]
    metricas_euler2 = Mt.run_metricas(Obs2, Sim2)
    # print 230 * '*'
    # print '///Metricas Euler con fpd02///'
    # print metricas_euler2.to_string()
    # print 230 * '*'

    # ***************************************** Metricas fdp01 vs fdp02 ***************************************** #
    comparativo = metricas_euler1.append(metricas_euler2, ignore_index=True)
    print 233 * '*'
    print '///Comparativo de metricas de desempeño///'
    print comparativo.to_string()
    print 233 * '*'

    # Ploteo de resultados
    # Galeria de estilo de ploteo
    plt.style.use('ggplot')
    # Linea a 45 grados
    L45 = [1, 790]

    # Subplot 1 - 1
    plt.subplot(221)
    datos['Q'].plot()
    plt.subplot(221)
    datos['Euler_fpd01'].plot()
    plt.legend(['Q', 'Euler fdp01'], loc='best')

    # Subplot 1 - 2
    plt.subplot(222)
    plt.plot(L45, L45, 'g')
    plt.scatter(datos['Q'], datos['Euler_fpd01'])
    plt.axis('equal')
    plt.title(r'Obs vs Sim. Euler fdp01 - ${r}^2$ = ' + str(np.corrcoef(datos['Q'], datos['Euler_fpd01'])[0, 1]))

    # Subplot 2 - 1
    plt.subplot(223)
    datos['Q'].plot()
    plt.subplot(223)
    datos['Euler_fpd02'].plot()
    plt.legend(['Q', 'Euler fdp02'], loc='best')

    # Subplot 2 - 2
    plt.subplot(224)
    plt.plot(L45, L45, 'g')
    plt.scatter(datos['Q'], datos['Euler_fpd02'])
    plt.title(r'Obs vs Sim. Euler fdp02 - ${r}^2$ = ' + str(np.corrcoef(datos['Q'], datos['Euler_fpd02'])[0, 1]))
    plt.axis('equal')
    # plt.tight_layout()
    plt.show()


