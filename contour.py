# -*- coding: utf-8 -*-
"""
|/------------Gamma Ray Burst Attenuation-----------\|
PHYS30302 -- Assignment 1 -- General Physics Poster
|/---------------------------------------------------------\|

The purpose of this document is to calculate the number of gamma ray photons
that are transmitted through the air to a known area at a known distance.


Bibliography:
    [1] - Wikipedia, "Nuclear Weapon Yield", https://en.wikipedia.org/wiki/Nuclear_weapon_yield
    [2] - Glasstone, S., and Dolan, P., J., "The Effects of Nuclear Weapons", 1977,  UNITED STATES
    DEPARTMENT OF DEFENSE and the ENERGY RESEARCH AND DEVELOPMENT ADMINISTRATION
    [3] - https://iba.online/knowledge/en/raeume-planen/office-planning/body-dimensions/ (Could defo get a better reference)

Created: Tue Mar 19 19:00:49 2024
Last Updated: Tue Mar 19 __:__:__ 2024

@author: Charlie Fynn Perkins, UID: 10839865 0

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Physical Constants
GAMMA_RAY_DIST = np.array([[0.75, 0.7], [2, 0.2], [4.5, 0.09], [8, 0.09], [
                          12, 0.02]])  # [Bin Max MeV, Proportion]

# Graph Related Globals
NUM_POINTS = 100
SAVE_DIR = "/Users/archiebrooks/Documents/Uni/Group-Project/genphys_poster/pictures"
FILENAME = "dose_distance_kt"

#approximating Harrison Ford to be a cuboid with measurements taken from [3]
HARRISON_THICKNESS = 0.35 #m
HARRISON_AREA = 0.925 #m^2
HARRISON_WEIGHT = 85 #kg

# Attenuation Values
AREA_TOP = 0.7186  #m^2
AREA_SIDE = 1.469 #m^2

DISTANCE_ARRAY = np.linspace(0,6.5,110)
YIELD_ARRAY = np.linspace(0,74,110)

# Bomb Setup
BOMB_YIELD = 60  # kt


def photon_yield(bomb_energy):
    """


    Parameters
    ----------
    bomb_energy : float, int
        Fission bomb TNT equivalent yield

    Returns
    -------
    photon_dist : np.array([int...])
        Number of photons at burst in gamma ray ranges in [2]

    """
    bomb_energy *= ((1/0.239) * 10**12)  # Convert to J from paragraph 1 in [1]
    instantaneous_gamma_energy = bomb_energy * 7/200  # From table 1.43 in [2]

    photon_dist = []
    for i in range(len(GAMMA_RAY_DIST)):
        # Get the median energy in the ith bin
        if i > 0:
            midpoint_energy = (
                GAMMA_RAY_DIST[i][0] - GAMMA_RAY_DIST[i-1][0]) / 2
            midpoint_energy += GAMMA_RAY_DIST[i][0]
        else:
            midpoint_energy = GAMMA_RAY_DIST[0][0] / 2

        # Calculate and append the number of photons from (bomb energy/photon energy)
        photon_dist.append(instantaneous_gamma_energy *
                           GAMMA_RAY_DIST[i][1] / (midpoint_energy * 1.60218e-13))
    return photon_dist


def photon_energy_dataset(photon_dist, total_points=1000):
    """


    Parameters
    ----------
    pphoton_dist : np.array([int...])
        Number of photons at burst in gamma ray ranges in [2]
    total_points : int, optional
        Number of points for the graph. The default is 1000.

    Returns
    -------
    energies : np.array([float])
        Numpy array of photon energies, MeV
    photons : np.array([float])
        Numpy array of number of photons corresponding to ith energy in energies

    """
    energies = np.array([])
    photons = np.array([])

    for i in range(len(photon_dist)):
        # Get the current and previous energies
        if i > 0:
            current_energy = GAMMA_RAY_DIST[i][0]
            previous_energy = GAMMA_RAY_DIST[i-1][0]
        else:
            current_energy = GAMMA_RAY_DIST[i][0]
            previous_energy = 0

        # Calculates number of points for space in graph
        point_proportion = (
            current_energy - previous_energy) / GAMMA_RAY_DIST[4][0]
        num_points = int(round(total_points * point_proportion, 0))

        # If this were Monte Carlo, the simulation could go here...
        new_energies = np.random.uniform(previous_energy,current_energy, num_points)
        energies = np.hstack((energies, new_energies))
        photons = np.hstack(
            (photons, np.full(num_points, photon_dist[i]/num_points)))
    return energies, photons


def attenuation_coeff_air(energy):
    return (8.36988619e-05) * np.power(energy, -4.60402706e-01)


def attenuation_coeff_lead_1(energy):
    return ((0.39830182) * np.power(energy, -1.63962241)) + 0.39926712

def attenuation_coeff_lead_2(energy):
    return 2.49309182e-04*(energy**3) + -6.71447586e-03*(energy**2) + 6.48613057e-02*energy + 3.24296337e-01

def attenuation_coeff_steel(energy):
    return 0.27869492 * np.power(energy, -0.81106709) + 0.17571693

def attenuation_coeff_water(energy):
    return 0.07424027*(np.power(energy,-0.46058569)) - 0.00462844

def attenuate(energies, photons, distance, material="air"):
    if material=="air":
        attenuated_photons = photons.tolist()
        for i in range(len(energies)):
            attenuated_photons[i] *= np.exp(-1 * attenuation_coeff_air(energies[i])
                                        * (distance * 100000))
    elif material=="steel":
        attenuated_photons = photons
        for i in range(len(energies)):
            attenuated_photons[i] *= np.exp(-1 * attenuation_coeff_steel(energies[i])
                                        * (distance * 100000))
    elif material=="lead":
        attenuated_photons = photons
        for i in range(len(energies)):
            if energies[i]<3:
                attenuated_photons[i] *= np.exp(-1 * attenuation_coeff_lead_1(energies[i])
                                        * (distance * 100000))
            else:
                attenuated_photons[i] *= np.exp(-1 * attenuation_coeff_lead_2(energies[i])
                                        * (distance * 100000))
    elif material=="water":
        attenuated_photons = photons
        for i in range(len(energies)):
            attenuated_photons[i] *= np.exp(-1 * attenuation_coeff_water(energies[i])
                                        * (distance * 100000))
    return (energies, attenuated_photons)


#this needs work
def graph(dose_array): 

    levels = [0.5, 1.5, 3.5, 8, 80]
    colors = ['orange', 'orangered', 'red', 'firebrick', 'maroon']
    labels = [f'Risk of Cancer', f'10% Fatality Rate', f'50% Fatality Rate', f'100% Fatality Rate', 'Instant Death']

    # Plot the contour map with contours and a logarithmic colorbar
    fig, ax = plt.subplots(figsize=(8, 8))
    contour = ax.contour(DISTANCE_ARRAY, YIELD_ARRAY, dose_array, levels=levels, colors=colors)
    
    ax.clabel(contour, fmt='%1.1f Sv', colors="k", inline=True, inline_spacing=20, fontsize=10)


    # Set axis labels and title
    ax.set_xlabel('Distance from explosion (km)')
    ax.set_ylabel('Yield of Explosion (kT)')
    ax.set_title('Contours of death')

    handles = [plt.Line2D([0], [0], linestyle='-', color=color, linewidth=2) for color in colors]
    ax.legend(handles, labels, loc='upper left')
    ax.set_ylim(-0.5)
    fig.tight_layout()
    fig.savefig(SAVE_DIR + "//" + FILENAME + ".png", dpi=800)
    fig.show()
    fig.clf()



if __name__ == "__main__":

    dose_array = np.zeros((len(DISTANCE_ARRAY),len(YIELD_ARRAY)))
    for i in range(len(YIELD_ARRAY)):
        energies, photons = photon_energy_dataset(
        photon_yield(YIELD_ARRAY[i]), NUM_POINTS)
        print(np.sum(energies*photons))
        for j in range(len(DISTANCE_ARRAY)):
            energies, air_photons = attenuate(energies, photons, DISTANCE_ARRAY[j])
            energies, steel_photons = attenuate(energies, air_photons, 10e-5, material="steel")
            energies, lead_photons = attenuate(energies, steel_photons, 3.35e-5, material="lead")

            incident_photons = np.sum(lead_photons)*HARRISON_AREA/(4*np.pi*(DISTANCE_ARRAY[j]+13.35e-5)**2)

            energies, water_photons = attenuate(energies, lead_photons, 0.35e-3, material="water")
            absorbed_photons = incident_photons-np.sum(water_photons)*HARRISON_AREA/(4*np.pi*(DISTANCE_ARRAY[j]+13.35e-5)**2)

            dose_array[i,j] = (np.sum(absorbed_photons*energies)*1.6e-13)/85

    graph(dose_array)