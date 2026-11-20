# %%
from four_vectors import FourVector
import data_parser as dp
import dal_raklev_equations as dr
import numpy as np
import random

valid_combs = [
    [dp.Hadron.PROTON, dp.Hadron.NEUTRON],
    [dp.Hadron.NEUTRON, dp.Hadron.PROTON],
]

DR_valid_combs = [
    [dp.Hadron.PROTON, dp.Hadron.NEUTRON],
    [dp.Hadron.NEUTRON, dp.Hadron.PROTON],
    [dp.Hadron.PROTON, dp.Hadron.PROTON],
    [dp.Hadron.NEUTRON, dp.Hadron.NEUTRON],
]

deut_mass = 1.8756129


def get_all_pairs(particles):
    pairs = []
    if pairs is not None:
        for i, p1 in enumerate(particles):
            for p2 in particles[i+1:]:
                combi = [p1.particle, p2.particle]
                if (p2.mat_type == p1.mat_type) and (combi in valid_combs):
                    pairs.append(dp.Pairing(p1.mat_type, np.array([p1, p2])))
        pairs_list = np.array(pairs)
        np.random.shuffle(pairs_list)
        return pairs_list
    else:
        return np.array([])


def find_k(particles):
    p1 = particles[0]
    p2 = particles[1]
    p_diff = (p1 - p2).p[1:]
    return np.sqrt(np.sum(np.power(p_diff, 2)))


def validate_pair(pair, cutoff):
    pair.com()
    k = find_k(pair.particles)
    pair.lab()
    if k <= cutoff:
        return True
    return False


def dr_validate(pair):
    pair.com()
    k = find_k(pair.particles)
    pair.lab()
    types = [pair.particles[0].particle, pair.particles[1].particle]
    sigmas = []
    sigmas.append(dr.pn_dg_cross(k))
    sigmas.append(dr.hh_pi0_cross(k, types))
    sigmas.append(dr.hh_dpi_cross(k))
    sigmas.append(dr.pn_pipm_cross(k))
    if (k < 0.3) or (k > 3):
        return 0
    else:
        sigmas = np.array(sigmas)
        weights = sigmas/np.sum(sigmas)
        total = np.sum(sigmas)
        sigmas /= total
        options = [1, 2, 3, 4]
        choice = random.choices(options, weights=weights)
        return choice


def get_valid_pairs(event, max_depth, cutoff):
    particles = event.particles
    pairs = get_all_pairs(particles)
    valid_pairs = []
    itr = 0
    while (pairs.size != 0) and (itr < max_depth):
        itr += 1
        for pair in pairs:
            if validate_pair(pair, cutoff):
                valid_pairs.append(pair)
                for part in valid_pairs[-1].particles:
                    particles = np.delete(particles, np.where(particles == part))
                pairs = get_all_pairs(particles)
                np.random.shuffle(pairs)
                break
    return np.array(valid_pairs)


def get_dr_valid_pairs(event, max_depth):
    particles = event.particles
    pairs = get_all_pairs(particles)
    valid_pairs = []
    itr = 0
    while (pairs.size != 0) and (itr < max_depth):
        itr += 1
        for pair in pairs:
            pathway = dr_validate(pair)
            if pathway != 0:
                pair.set_pathway(pathway)
                valid_pairs.append(pair)
                for part in valid_pairs[-1].particles:
                    particles = np.delete(particles, np.where(particles == part))
                pairs = get_all_pairs(particles)
                np.random.shuffle(pairs)
                break
    return np.array(valid_pairs)


def form_dr_deuterons(event, max_depth=50):
    valid_pairs = get_dr_valid_pairs(event, max_depth)
    deuterons = []
    if valid_pairs.size != 0:
        for pair in valid_pairs:
            pair.com()
            X = pair.get_inter_particle()
            pair.lab()
            invar_mass = X.p[0]
            p_mag = ((invar_mass**2)-(deut_mass**2))/(2*invar_mass)
            p_vect = p_mag * isotropic_sample()
            E = np.sqrt((invar_mass**2) + (np.sum(p_vect)**2))
            fv = FourVector([E, *p_vect])
            d = dp.Deuteron(pair.mat_type, fv, pair.cm_velo, pair.pathway)
            d.lab()
            deuterons.append(d)
        return deuterons
    return None


def isotropic_sample():
    """
        Returns a uniformly random direction on a sphere

        Returns
        -------
        vct : `array`
            the coordinates of the unit vector (i.e. from origin to `vct`) of form `[x,y,z]`
    """
    phi = np.random.uniform(0, 2*np.pi)
    theta = np.arccos(np.random.uniform(-1, 1))
    coords = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
    return coords


def form_deuterons(event, max_depth=50, cutoff=0.1):
    valid_pairs = get_valid_pairs(event, max_depth, cutoff)
    deuterons = []
    if valid_pairs.size != 0:
        for pair in valid_pairs:
            pair.com()
            X = pair.get_inter_particle()
            pair.lab()
            invar_mass = X.p[0]
            p_mag = ((invar_mass**2)-(deut_mass**2))/(2*invar_mass)
            p_vect = p_mag * isotropic_sample()
            E = np.sqrt((invar_mass**2) + (np.sum(p_vect)**2))
            fv = FourVector([E, *p_vect])
            d = dp.Deuteron(pair.mat_type, fv, pair.cm_velo)
            d.lab()
            deuterons.append(d)
        return deuterons
    return None

# %%
