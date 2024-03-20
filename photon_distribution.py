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

# Physical Constants
GAMMA_RAY_DIST = np.array([[0.75, 0.7], [2, 0.2], [4.5, 0.09], [8, 0.09], [
                          12, 0.02]])  # [Bin Max MeV, Proportion]

# Graph Related Globals
NUM_POINTS = 100
SAVE_DIR = "/Users/archiebrooks/Documents/Uni/Group-Project/genphys_poster/pictures"
FILENAME = "energy_distribution"

# Attenuation Values
AREA_TOP = 0.7186  #m^2
AREA_SIDE = 1.469 #m^2
DISTANCE = 5  # km

# Bomb Setup
BOMB_YIELD = 100  # kt


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

        new_energies = np.linspace(previous_energy,current_energy, num_points)
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
    else:
        print("Material not recognised")
        exit
    return (energies, attenuated_photons)


def graph(energies, photons, attenuated_photons):
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))

    ax.set(title="Attenuated Gamma-Ray Photons from " +
           str(BOMB_YIELD) + " kt Nuclear Blast")
    ax.set(xlabel="Gamma Ray Energy, MeV")
    ax.set(ylabel="Number of Photons")

    ax.plot(energies, photons, linestyle="dashed", label="At Burst")
    ax.plot(energies, attenuated_photons, linestyle="solid", label="At " + str(DISTANCE) + " km")

    ax.set_yscale("log")
    ax.legend()
    ax.grid(alpha=0.5)
    ax.set_ylim(1e-2,1e26)

    fig.tight_layout()
    fig.savefig(SAVE_DIR + "//" + FILENAME + ".png", dpi=800)
    fig.show()
    fig.clf()


if __name__ == "__main__":
    energies, photons = photon_energy_dataset(
        photon_yield(BOMB_YIELD), NUM_POINTS)

    energies, air_photons = attenuate(energies, photons, DISTANCE)

    # Here goes code to reduce number of photons in the solid angle of the fridge

    graph(energies, photons, air_photons)

    print("Total number of photons at distance", DISTANCE, "km")
    print(np.sum(air_photons))

    print("As percentage of initial...")
    print(np.sum(air_photons)/np.sum(photons) * 100, "%")

    print("Total number of photons incident on fridge at ", DISTANCE, " km:")
    print(np.sum(air_photons)*AREA_SIDE/(4*np.pi*DISTANCE**2))

    energies, steel_photons = attenuate(energies, air_photons, 10e-5, material="steel")
    energies, lead_photons = attenuate(energies, steel_photons, 3.35e-5, material="lead")

    print("Total number of photons which pass through fridge wall at ", DISTANCE, " km:")
    print(np.sum(lead_photons)*AREA_SIDE/(4*np.pi*DISTANCE**2))

    print("As a percentage:")
    print((np.sum(lead_photons)*AREA_SIDE/(4*np.pi*DISTANCE**2))/np.sum(photons) * 100)


    FILENAME = "energy_after_fridge"
    graph(energies, photons, lead_photons)
    
