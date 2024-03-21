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

Created: Tue Mar 19 19:00:49 2024
Last Updated: Tue Mar 19 __:__:__ 2024

@author: Charlie Fynn Perkins, UID: 10839865 0
"""

import numpy as np
import matplotlib.pyplot as plt
import random

# Physical Constants
GAMMA_RAY_DIST = np.array([[0.75, 0.7], [2, 0.2], [4.5, 0.09], [8, 0.09], [
                          12, 0.02]])  # [Bin Max MeV, Proportion]

# Graph Related Globals
NUM_POINTS = 1000
SAVE_DIR = "/Users/archiebrooks/Documents/Uni/Group-Project/genphys_poster/pictures"
FILENAME = "fridge_transport"

# Attenuation Values
AREA_TOP = 0.7186  #m^2
AREA_SIDE = 1.469 #m^2
DISTANCE = 5 # km

DISTANCE_ARRAY = np.linspace(1,15,500)

# Bomb Setup
BOMB_YIELD = 10  # kt


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


def graph(distances, air, steel, lead):
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))

    ax.set(title="Gamma ray transport through fridge walls at various distances from a " +
           str(BOMB_YIELD) + " kt Nuclear Blast")
    ax.set(xlabel="Distance from explosion (km)")
    ax.set(ylabel="Number of Photons")

    ax.plot(distances, air, linestyle="solid",
            color="indigo", label="Incident")
    ax.plot(distances, steel, linestyle="solid",
            color="darkgoldenrod", label="After steel layer")
    ax.plot(distances, lead, linestyle="solid",
            color="magenta", label="After lead layer")
    ax.set_yscale("log")
    ax.set_ylim(10e-3)
    ax.legend()
    ax.grid(alpha=0.5)

    fig.tight_layout()
    fig.savefig(SAVE_DIR + "//" + FILENAME + ".png", dpi=800)
    fig.show()
    fig.clf()


if __name__ == "__main__":
    energies, photons = photon_energy_dataset(
        photon_yield(BOMB_YIELD), NUM_POINTS)

    photon_array = []
    air_array = []
    steel_array = []
    lead_array = []
    
    for d in DISTANCE_ARRAY:
        air_energies, air_photons = attenuate(energies, photons, d)
        air_results = np.sum(air_photons)*AREA_SIDE/(4*np.pi*d**2)

        steel_energies, steel_photons = attenuate(air_energies, air_photons, 10e-5, material="steel")
        steel_results = np.sum(steel_photons)*AREA_SIDE/(4*np.pi*(d+10e-5)**2)

        lead_energies, lead_photons = attenuate(steel_energies, steel_photons, 3.35e-5, material="lead")
        lead_results = np.sum(lead_photons)*AREA_SIDE/(4*np.pi*(d+13.35e-5)**2)

        air_array.append(air_results)
        steel_array.append(steel_results)
        lead_array.append(lead_results)
    
    graph(DISTANCE_ARRAY, air_array, steel_array, lead_array)