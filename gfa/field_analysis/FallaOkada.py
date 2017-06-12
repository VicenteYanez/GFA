#! /usr/bin/env python

"""
@autor: Klaus Bataille
@editor: vicente
Clase para manejar fallas de Okada
Input:
param[
      D(profundidad a punto origen de okada),
      W(longitud de falla desde el origen hacia superficie),
      largo,
      rake,
      strike
      dip
      y_origen
      x_origen
      ]
malla[x,y]
"""

from math import *
from numpy import *
import scipy as scp
from scipy import optimize
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata


class FallaOkada():
    def __init__(self, parametros, malla):
        # desgloce parametros
        self.D = parametros[0]
        self.W = parametros[1]
        self.largo = parametros[2]
        self.rake = parametros[3]
        self.strike = parametros[4]
        self.dip = parametros[5]
        y_origen = parametros[6]
        x_origen = parametros[7]

        # determina el origen del sistema de Okada
        lat0, lon0, alpha21 = self._vinc_pt(y_origen, x_origen,
                                           degrees(self.strike+pi/2.),
                                           self.W*cos(self.dip))
        # conversion malla a coordenadas de okada
        self.x, self.y = self._proj_mesh2okada(malla[1], malla[0], lat0, lon0)

    def tensor_deformacion(self):
        p = self.y*cos(self.dip) + self.D*sin(self.dip)
        q = self.y*sin(self.dip) - self.D*cos(self.dip)
        e = array([self.x, self.x, self.x-self.largo, self.x-self.largo]).T
        eta = array([p, p-self.W, p, p-self.W]).T
        qq = array([q, q, q, q]).T    # b = 4

        ytg = eta*cos(self.dip) + qq*sin(self.dip)
        dtg = eta*sin(self.dip) - qq*cos(self.dip)
        R = power(e**2 + eta**2 + qq**2, 0.5)
        X = power(e**2 + qq**2, 0.5)

        # multiplicacion *0.5 equivale a u/(lambda+u)
        # condicionales que dependen del manteo de la falla
        if degrees(self.dip) != 90:
            K1 = 0.5*e/cos(self.dip)*(1./R/(R+dtg) - sin(self.dip)/R/(R+eta))
            K3 = 0.5/cos(self.dip)*(qq/R/(R+eta) - ytg/R/(R+dtg))

            J1 = (0.5/cos(self.dip)*(e*e/R/(R+dtg)**2 - 1/(R+dtg)) -
                 sin(self.dip)*K3/cos(self.dip))
            J2 = (0.5/cos(self.dip)*(e*ytg/R/(R+dtg)**2) -
                 sin(self.dip)*K1/cos(self.dip))

        elif degrees(self.dip) == 90:
            K1 = 0.5*e*qq/(R*(R+dtg)**2)
            K3 = 0.5*sin(self.dip)/(R+dtg)*(e**2/(R*(R+dtg))-1.)

            J1 = 0.5*qq/(2*(R+dtg)**2)*(2*e**2/(R*(R+dtg))-1.)
            J2 = 0.5*e*sin(self.dip)/(2*(R+dtg)**2)*(2*qq**2/(R*(R+dtg)) - 1.)

        K4 = (2.*ytg/R/(R+e) + e*cos(self.dip)/R/(R+eta)) * sin(self.dip)
        K2 = 0.5*(- sin(self.dip) / R + qq*cos(self.dip)/R/(R+eta)) - K3
        J3 = - 0.5*e/R/(R+eta) - J2
        J4 = -0.5*(cos(self.dip)/R + qq*sin(self.dip)/R/(R+eta)) - J1

        Apsi = (2*R+e)/(R**3)/((R+e)**2)
        Aeta = (2*R+eta)/(R**3)/((R+eta)**2)

        # self.dip-slip
        ux_dx_ds = sin(self.rake)/(2*pi) * (e*qq/(R**3) + J3*sin(self.dip)*cos(self.dip))
        ux_dy_ds = sin(self.rake)/(2*pi) * (ytg*qq/(R**3) - sin(self.dip)/R +
                                       J1*sin(self.dip)*cos(self.dip))
        uy_dx_ds = sin(self.rake)/(2*pi) * (ytg*qq/(R**3) + qq*cos(self.dip)/R/(R+eta) +
                                       J1*sin(self.dip)*cos(self.dip))
        uy_dy_ds = sin(self.rake)/(2*pi) * (ytg*ytg*qq*Apsi - K4 + J2*sin(self.dip)*cos(self.dip))

        uz_dx_ds = sin(self.rake)/(2*pi) * (dtg*qq/(R**3) + qq*sin(self.dip)/R/(R+eta) +
                                       K3*sin(self.dip)*cos(self.dip))

        K4 = (2*dtg/R/(R+e) + e*sin(self.dip)/R/(R+eta)) * sin(self.dip)
        uz_dy_ds = sin(self.rake)/(2*pi) * (ytg*dtg*qq*Apsi - K4 + K1*sin(self.dip)*cos(self.dip))

        # strike-slip
        ux_dx_ss = cos(self.rake)/(2*pi) * (e*e*qq*Aeta - J1*sin(self.dip))
        ux_dy_ss = cos(self.rake)/(2*pi) * (e*e*e*dtg/(R**3)/(eta*eta+qq*qq) -
                                       sin(self.dip)*(e*e*e*Aeta+J2))
        uy_dx_ss = cos(self.rake)/(2*pi) * (e*qq*cos(self.dip)/(R**3) +
                                       sin(self.dip)*(e*qq*qq*Aeta - J2))

        K4 = (qq**3)*Aeta*sin(self.dip) - 2*qq*sin(self.dip)/R/(R+eta) - (e*e+eta*eta)*cos(self.dip)/(R**3) - J4
        K4 = K4 * sin(self.dip) + ytg*qq*cos(self.dip)/(R**3)
        uy_dy_ss = cos(self.rake)/(2*pi) * K4

        uz_dx_ss = cos(self.rake)/(2*pi) * (-e*qq*qq*Aeta*cos(self.dip) + (e*qq/(R**3) -
                                       K1)*sin(self.dip))
        K5 = dtg*qq/(R**3)*cos(self.dip) + (e*e*qq*Aeta*cos(self.dip)-sin(self.dip)/R +
                                       ytg*qq/(R**3)-K2)*sin(self.dip)
        uz_dy_ss = cos(self.rake)/(2*pi) * K5

        # representacion chinnery self.dip-slip
        ux_dx_d = ux_dx_ds.T[0]-ux_dx_ds.T[1]-ux_dx_ds.T[2]+ux_dx_ds.T[3]
        ux_dy_d = ux_dy_ds.T[0]-ux_dy_ds.T[1]-ux_dy_ds.T[2]+ux_dy_ds.T[3]
        uy_dx_d = uy_dx_ds.T[0]-uy_dx_ds.T[1]-uy_dx_ds.T[2]+uy_dx_ds.T[3]
        uy_dy_d = uy_dy_ds.T[0]-uy_dy_ds.T[1]-uy_dy_ds.T[2]+uy_dy_ds.T[3]
        uz_dx_d = uz_dx_ds.T[0]-uz_dx_ds.T[1]-uz_dx_ds.T[2]+uz_dx_ds.T[3]
        uz_dy_d = uz_dy_ds.T[0]-uz_dy_ds.T[1]-uz_dy_ds.T[2]+uz_dy_ds.T[3]
        # uzd = uz_ds.T[0]-uz_ds.T[1]-uz_ds.T[2]+uz_ds.T[3]

        # representacion chinnery strike-slip
        ux_dx_s = ux_dx_ss.T[0]-ux_dx_ss.T[1]-ux_dx_ss.T[2]+ux_dx_ss.T[3]
        ux_dy_s = ux_dy_ss.T[0]-ux_dy_ss.T[1]-ux_dy_ss.T[2]+ux_dy_ss.T[3]
        uy_dx_s = uy_dx_ss.T[0]-uy_dx_ss.T[1]-uy_dx_ss.T[2]+uy_dx_ss.T[3]
        uy_dy_s = uy_dy_ss.T[0]-uy_dy_ss.T[1]-uy_dy_ss.T[2]+uy_dy_ss.T[3]
        uz_dx_s = uz_dx_ss.T[0]-uz_dx_ss.T[1]-uz_dx_ss.T[2]+uz_dx_ss.T[3]
        uz_dy_s = uz_dy_ss.T[0]-uz_dy_ss.T[1]-uz_dy_ss.T[2]+uz_dy_ss.T[3]

        # solucion
        ux_dx = ux_dx_d + ux_dx_s
        ux_dy = ux_dy_d + ux_dy_s
        uy_dx = uy_dx_d + uy_dx_s
        uy_dy = uy_dy_d + uy_dy_s
        uz_dx = uz_dx_d + uz_dx_s
        uz_dy = uz_dy_d + uz_dy_s

        # at free surface:
        ux_dz = - uz_dx
        uy_dz = - uz_dy
        uz_dz = -1./3. * (ux_dx + uy_dy)

        return ux_dx, ux_dy, uy_dx, uy_dy, uz_dx, uz_dy, uz_dz


    # OUTPUT: ux,uy,uz en coordenadas de okada
    def desplaz(self):
        # notacion  de Chinnery:f(e,eta)||= f(x,p)-f(x,p-W)-f(x-L,p)+f(x-L,W-p)
        p = self.y*cos(self.dip) + self.D*sin(self.dip)
        q = self.y*sin(self.dip) - self.D*cos(self.dip)
        e = array([self.x, self.x, self.x - self.largo, self.x - self.largo]).T
        eta = array([p, p-self.W, p, p-self.W]).T
        qq = array([q, q, q, q]).T    # b = 4

        ytg = eta*cos(self.dip) + qq*sin(self.dip)
        dtg = eta*sin(self.dip) - qq*cos(self.dip)
        R = power(e**2 + eta**2 + qq**2, 0.5)
        X = power(e**2 + qq**2, 0.5)

        if degrees(self.dip) != 90:
            I5 = (1/cos(self.dip))*scp.arctan((eta*(X+qq*cos(self.dip)) +
                                              X*(R+X)*sin(self.dip))/
                                              (e*(R+X)*cos(self.dip)))

            I4 = .5/cos(self.dip)*(scp.log(R+dtg)-sin(self.dip)*scp.log(R+eta))

            I1 = (.5*((-1./cos(self.dip))*(e/(R+dtg))) -
                 (sin(self.dip)*I5/cos(self.dip)))

            I3 = (.5*(1/cos(self.dip)*(ytg/(R+(dtg))) -
                 scp.log(R+eta))+(sin(self.dip) * I4/cos(self.dip)))

        if degrees(self.dip) == 90:
            I5 = -.5*e*sin(self.dip)/(R+dtg)
            I4 = -.5*qq/(R+dtg)
            I3 = .25*(eta/(R+dtg) + ytg/(R+dtg)**2 - scp.log(R+eta))
            I1 = -.25*e*qq/(R+dtg)**2

        I2 = 0.5*(-scp.log(R+eta))-I3

        # self.dip-slip
        ux_ds = -sin(self.rake)/(2*pi)*(qq/R-I3*sin(self.dip)*cos(self.dip))
        uy_ds = -sin(self.rake)/(2*pi)*((ytg*qq/R/(R+e)) +
                                   (cos(self.dip)*scp.arctan(e*eta/qq/R)) -
                                   (I1*sin(self.dip)*cos(self.dip)))
        uz_ds = -sin(self.rake)/(2*pi)*((dtg*qq/R/(R+e)) +
                                   (sin(self.dip)*scp.arctan(e*eta/qq/R)) -
                                   (I5*sin(self.dip)*cos(self.dip)))

        # strike-slipe
        ux_ss = -cos(self.rake)/(2*pi)*((e*qq/R/(R+eta)) +
                                   (scp.arctan(e*eta/(qq*R)))+I1*sin(self.dip))
        uy_ss = -cos(self.rake)/(2*pi)*((ytg*qq/R/(R+eta)) +
                                   qq*cos(self.dip)/(R+eta)+I2*sin(self.dip))
        uz_ss = -cos(self.rake)/(2*pi)*((dtg*qq/R/(R+eta)) +
                                   qq*sin(self.dip)/(R+eta)+I4*sin(self.dip))

        # representacion chinnery self.dip-slip
        uxd = ux_ds.T[0]-ux_ds.T[1]-ux_ds.T[2]+ux_ds.T[3]
        uyd = uy_ds.T[0]-uy_ds.T[1]-uy_ds.T[2]+uy_ds.T[3]
        uzd = uz_ds.T[0]-uz_ds.T[1]-uz_ds.T[2]+uz_ds.T[3]

        # representacion chinnery strike-slip
        uxs = ux_ss.T[0]-ux_ss.T[1]-ux_ss.T[2]+ux_ss.T[3]
        uys = uy_ss.T[0]-uy_ss.T[1]-uy_ss.T[2]+uy_ss.T[3]
        uzs = uz_ss.T[0]-uz_ss.T[1]-uz_ss.T[2]+uz_ss.T[3]

        # cantidad de desplazamiento
        uxs = uxs
        uys = uys
        uzs = uzs

        uxd = uxd
        uyd = uyd
        uzd = uzd

        # suma componentes strike y dip slip.
        ux = uxd + uxs
        uy = uyd + uys
        uz = uzd + uzs

        # proyeccion a las componentes geograficas
        Ue = ux*sin(self.strike) - uy*cos(self.strike)
        Un = ux*cos(self.strike) + uy*sin(self.strike)

        # para revisar valores
        if False:
            print(ux, uy, uz)

        return Ue, Un, uz


###############################################################################
# OUTPUT: coordenadas xi,yi proyectadas en los ejes de Okada
    def _proj_mesh2okada(self, lat, lon, lat0, lon0):
        if len(lat) and len(lon) > 1:
            # lat lon, rango en la falla y
            # lat0 lon0 sistema de origen de los ejes de okada
            dlat = zeros(len(lat))
            dlon = zeros(len(lon))
            for i in range(len(lat)):
                llat = lat[i]
                llon = lon[i]
                l, az, baz = self._vinc_dist(lat0, lon0, llat, llon)
                dlat[i] = l*cos(radians(az))
                dlon[i] = l*sin(radians(az))

            Dla = array(dlat)
            Dlo = array(dlon)

            xi = Dla*cos(self.strike) + Dlo*sin(self.strike)
            yi = Dla*sin(self.strike) - Dlo*cos(self.strike)
        else:
            print("Malla de un solo punto")
            xi = lon
            yi = lat

        return xi, yi


###############################################################################
# Transforma ux,uy (coordenadas okada) a Ue,Un (coordenadas espaciales)
    def _okada2cardinals(self, ux, uy):
        Ue = ux*sin(self.strike) - uy*cos(self.strike)
        Un = ux*cos(self.strike) + uy*sin(self.strike)
        return Ue, Un

#
# ---------------------------------------------------------------------
# |                                                                     |
# |	geodetic.cc -  a collection of geodetic functions                   |
# |	Jim Leven  - Dec 99                                                 |
# |                                                                     |
# | originally from:                                                    |
# | http://wegener.mechanik.tu-darmstadt.de/GMT-Help/Archiv/att-8710/Geodetic_py |                                                                   |
# |                                                                     |
# ---------------------------------------------------------------------
#
#
# ----------------------------------------------------------------------
# | Algrothims from Geocentric Datum of Australia Technical Manual	    |
# | 								                                    |
# | http://www.anzlic.org.au/icsm/gdatum/chapter4.html	        		|
# | 								                                    |
# | This page last updated 11 May 1999 	                				|
# | 								                                    |
# | Computations on the Ellipsoid	                    				|
# | 								                                    |
# | There are a number of formulae that are available           		|
# | to calculate accurate geodetic positions, 		            		|
# | azimuths and distances on the ellipsoid.			                |
# | 								                                    |
# | Vincenty's formulae (Vincenty, 1975) may be used 		            |
# | for lines ranging from a few cm to nearly 20,000 km, 	            |
# | with millimetre accuracy. 					                        |
# | The formulae have been extensively tested 		                    |
# | for the Australian region, by comparison with results       		|
# | from other formulae (Rainsford, 1955 & Sodano, 1965). 	            |
# |								                                        |
# | * Inverse problem: azimuth and distance from known 	        		|
# |			latitudes and longitudes 			                        |
# | * Direct problem: Latitude and longitude from known 	            |
# |			position, azimuth and distance. 		                    |
# | * Sample data 						                                |
# | * Excel spreadsheet 			                            		|
# | 								                                    |
# | Vincenty's Inverse formulae				                    		|
# | Given: latitude and longitude of two points                 		|
# |			(phi1, lembda1 and phi2, lembda2), 	                        |
# | Calculate: the ellipsoidal distance (s) and 	            		|
# | forward and reverse azimuths between the points (alpha12, alpha21).	|
# |									                                    |
# ----------------------------------------------------------------------

    def _vinc_dist(self, phi1,  lembda1,  phi2,  lembda2):
        """
        Returns the distance between two geographic points on the ellipsoid
        and the forward and reverse azimuths between these points.
        lats, longs and azimuths are in decimal degrees, distance in metres

        Returns ( s, alpha12,  alpha21 ) as a tuple
        """

        f = 1.0 / 298.257223563		# WGS84
        a = 6378137.0 			# metres
        if (abs(phi2 - phi1) < 1e-8) and (abs(lembda2 - lembda1) < 1e-8):
            return 0.0, 0.0, 0.0

        piD4 = math.atan(1.0)
        two_pi = piD4 * 8.0

        phi1 = phi1 * piD4 / 45.0
        lembda1 = lembda1 * piD4 / 45.0		# unfortunately lambda is a key word!
        phi2 = phi2 * piD4 / 45.0
        lembda2 = lembda2 * piD4 / 45.0

        b = a * (1.0 - f)

        TanU1 = (1-f) * tan(phi1)
        TanU2 = (1-f) * tan(phi2)

        U1 = atan(TanU1)
        U2 = atan(TanU2)

        lembda = lembda2 - lembda1
        last_lembda = -4000000.0		# an impossibe value
        omega = lembda

        # Iterate the following equations,
        #  until there is no significant change in lembda

        while (last_lembda < -3000000.0 or lembda != 0 and abs((last_lembda - lembda)/lembda) > 1.0e-9):

            sqr_sin_sigma = pow(cos(U2) * sin(lembda), 2) + pow((cos(U1) * sin(U2) - \
                    sin(U1) *  cos(U2) * cos(lembda) ), 2 )

            Sin_sigma = sqrt( sqr_sin_sigma )

            Cos_sigma = sin(U1) * sin(U2) + cos(U1) * cos(U2) * cos(lembda)

            sigma = atan2( Sin_sigma, Cos_sigma )

            Sin_alpha = cos(U1) * cos(U2) * sin(lembda) / sin(sigma)
            alpha = asin( Sin_alpha )

            Cos2sigma_m = cos(sigma) - (2 * sin(U1) * sin(U2) / pow(cos(alpha), 2) )

            C = (f/16) * pow(cos(alpha), 2) * (4 + f * (4 - 3 * pow(cos(alpha), 2)))

            last_lembda = lembda

            lembda = omega + (1-C) * f * sin(alpha) * (sigma + C * sin(sigma) * \
                    (Cos2sigma_m + C * cos(sigma) * (-1 + 2 * pow(Cos2sigma_m, 2) )))

        u2 = pow(cos(alpha),2) * (a*a-b*b) / (b*b)

        A = 1 + (u2/16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))

        B = (u2/1024) * (256 + u2 * (-128+ u2 * (74 - 47 * u2)))

        delta_sigma = B * Sin_sigma * (Cos2sigma_m + (B/4) * \
                (Cos_sigma * (-1 + 2 * pow(Cos2sigma_m, 2) ) - \
                (B/6) * Cos2sigma_m * (-3 + 4 * sqr_sin_sigma) * \
                (-3 + 4 * pow(Cos2sigma_m,2 ) )))

        s = b * A * (sigma - delta_sigma)

        alpha12 = atan2( (cos(U2) * sin(lembda)), \
                (cos(U1) * sin(U2) - sin(U1) * cos(U2) * cos(lembda)))

        alpha21 = atan2( (cos(U1) * sin(lembda)), \
                (-sin(U1) * cos(U2) + cos(U1) * sin(U2) * cos(lembda)))

        if ( alpha12 < 0.0 ) :
                alpha12 =  alpha12 + two_pi
        if ( alpha12 > two_pi ) :
                alpha12 = alpha12 - two_pi

        alpha21 = alpha21 + two_pi / 2.0
        if ( alpha21 < 0.0 ) :
                alpha21 = alpha21 + two_pi
        if ( alpha21 > two_pi ) :
                alpha21 = alpha21 - two_pi

        alpha12    = alpha12    * 45.0 / piD4
        alpha21    = alpha21    * 45.0 / piD4
        return s, alpha12,  alpha21

# END of Vincenty's Inverse formulae


# -----------------------------------------------------------------------------
# Vincenty's Direct formulae							|
# Given: latitude and longitude of a point (phi1, lembda1) and 			|
# the geodetic azimuth (alpha12) 						|
# and ellipsoidal distance in metres (s) to a second point,			|
# 										|
# Calculate: the latitude and longitude of the second point (phi2, lembda2) |
# and the reverse azimuth (alpha21).						|
# 										|
# -----------------------------------------------------------------------------

    def _vinc_pt(self, phi1, lembda1, alpha12, s):
        """
        Returns the lat and long of projected point and reverse azimuth
        given a reference point and a distance and azimuth to project.
        lats, longs and azimuths are passed in decimal degrees

        Returns ( phi2,  lambda2,  alpha21 ) as a tuple

        """

        f = 1.0 / 298.257223563		# WGS84
        a = 6378137.0 			# metres
        piD4 = atan(1.0)
        two_pi = piD4 * 8.0

        phi1 = phi1 * piD4 / 45.0
        lembda1 = lembda1 * piD4 / 45.0
        alpha12 = alpha12 * piD4 / 45.0
        if (alpha12 < 0.0):
                alpha12 = alpha12 + two_pi
        if ( alpha12 > two_pi ) :
                alpha12 = alpha12 - two_pi

        b = a * (1.0 - f)

        TanU1 = (1-f) * tan(phi1)
        U1 = atan( TanU1 )
        sigma1 = atan2( TanU1, cos(alpha12) )
        Sinalpha = cos(U1) * sin(alpha12)
        cosalpha_sq = 1.0 - Sinalpha * Sinalpha

        u2 = cosalpha_sq * (a * a - b * b ) / (b * b)
        A = 1.0 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * \
                (320 - 175 * u2) ) )
        B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2) ) )

        # Starting with the approximation
        sigma = (s / (b * A))

        last_sigma = 2.0 * sigma + 2.0	# something impossible

        # Iterate the following three equations
        #  until there is no significant change in sigma

        # two_sigma_m , delta_sigma
        while (abs((last_sigma - sigma) / sigma) > 1.0e-9):
                two_sigma_m = 2 * sigma1 + sigma

                delta_sigma = B * sin(sigma) * (cos(two_sigma_m) +
                                                (B/4) * (cos(sigma) *
                        (-1 + 2 * pow( cos(two_sigma_m), 2) -
                        (B/6) * cos(two_sigma_m) *
                        (-3 + 4 * pow(sin(sigma), 2 )) *
                        (-3 + 4 * pow( cos (two_sigma_m), 2 )))))

                last_sigma = sigma
                sigma = (s / (b * A)) + delta_sigma

        phi2 = atan2((sin(U1) * cos(sigma) + cos(U1) * sin(sigma) *
                     cos(alpha12)),
                     ((1-f) * sqrt(pow(Sinalpha, 2) + pow(sin(U1) *
                      sin(sigma) - cos(U1) * cos(sigma) * cos(alpha12), 2))))

        lembda = atan2( (sin(sigma) * sin(alpha12 )), (cos(U1) * cos(sigma) -  \
                sin(U1) *  sin(sigma) * cos(alpha12)))

        C = (f/16) * cosalpha_sq * (4 + f * (4 - 3 * cosalpha_sq ))

        omega = lembda - (1-C) * f * Sinalpha *  \
                (sigma + C * sin(sigma) * (cos(two_sigma_m) + \
                C * cos(sigma) * (-1 + 2 * pow(cos(two_sigma_m),2) )))

        lembda2 = lembda1 + omega

        alpha21 = atan2 ( Sinalpha, (-sin(U1) * sin(sigma) +  \
                cos(U1) * cos(sigma) * cos(alpha12)))

        alpha21 = alpha21 + two_pi / 2.0
        if ( alpha21 < 0.0 ) :
                alpha21 = alpha21 + two_pi
        if ( alpha21 > two_pi ) :
                alpha21 = alpha21 - two_pi

        phi2       = phi2       * 45.0 / piD4
        lembda2    = lembda2    * 45.0 / piD4
        alpha21    = alpha21    * 45.0 / piD4

        return phi2,  lembda2,  alpha21

  # END of Vincenty's Direct formulae

#--------------------------------------------------------------------------
# Notes:
#
# * "The inverse formulae may give no solution over a line
# 	between two nearly antipodal points. This will occur when
# 	lembda ... is greater than pi in absolute value". (Vincenty, 1975)
#
# * In Vincenty (1975) L is used for the difference in longitude,
# 	however for consistency with other formulae in this Manual,
# 	omega is used here.
#
# * Variables specific to Vincenty's formulae are shown below,
# 	others common throughout the manual are shown in the Glossary.
#
#
# alpha = Azimuth of the geodesic at the equator
# U = Reduced latitude
# lembda = Difference in longitude on an auxiliary sphere (lembda1 & lembda2
# 		are the geodetic longitudes of points 1 & 2)
# sigma = Angular distance on a sphere, from point 1 to point 2
# sigma1 = Angular distance on a sphere, from the equator to point 1
# sigma2 = Angular distance on a sphere, from the equator to point 2
# sigma_m = Angular distance on a sphere, from the equator to the
# 		midpoint of the line from point 1 to point 2
# u, A, B, C = Internal variables
#
#

# ******************************************************************
