# %%
import data_parser as dp
import numpy as np

valid_combs = [
    [dp.Nucleon.PROTON, dp.Nucleon.NEUTRON],
    [dp.Nucleon.NEUTRON, dp.Nucleon.PROTON],
]

DR_valid_combs = [
    [dp.Nucleon.PROTON, dp.Nucleon.NEUTRON],
    [dp.Nucleon.NEUTRON, dp.Nucleon.PROTON],
    [dp.Nucleon.PROTON, dp.Nucleon.PROTON],
    [dp.Nucleon.NEUTRON, dp.Nucleon.NEUTRON],
]


def get_all_pairs(particles):
    pairs = []
    if pairs is not None:
        for i, p1 in enumerate(particles):
            for p2 in particles[i+1:]:
                combi = [p1.particle, p2.particle]
                if (p2.mat_type == p1.mat_type) and (combi in valid_combs):
                    pairs.append(dp.Pairing(p1.mat_type, np.array([p1, p2])))
        return np.array(pairs)
    else:
        return np.array([])


def find_k(particles):
    p1 = particles[0]
    p2 = particles[1]
    p_diff = (p1 - p2).p[1:]
    return np.sqrt(np.sum(np.power(p_diff, 2)))


def validate_pair(pair, cutoff=0.2):
    pair.shift()
    k = find_k(pair.particles)
    if k >= cutoff:
        return True
    return False


def get_valid_pairs(event):
    particles = event.particles
    pairs = get_all_pairs(particles)
    valid_pairs = []
    while pairs.size != 0:
        for pair in pairs:
            if validate_pair(pair):
                valid_pairs.append(pair)
                for part in valid_pairs[-1].particles:
                    # print(part)
                    # print(valid_pairs[-1])
                    # print(valid_pairs[-1].particles)
                    # print(np.where(particles == part))
                    particles = np.delete(particles, np.where(particles == part))
                pairs = get_all_pairs(particles)
                break
    return np.array(valid_pairs)

# %%
