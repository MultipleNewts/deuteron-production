# %%
from dataclasses import dataclass, field
from enum import Enum
import numpy as np
from four_vectors import FourVector


# define a matter type class
class MatterType(Enum):
    MATTER = "m"
    ANTIMATTER = "a"


# define a nucleon type class
class Hadron(Enum):
    PROTON = "p"
    NEUTRON = "n"
    DEUTERON = "D"


# A particle dataclass for the datastructure
# the particle string starts with either "m-" or "a-" for matter and antimatter respectively
# then it is succeeded by "-p" or "-n" for proton or neutron respectively
# it also contains the corresponds four-vector and mass
# Note: energy and mass-energy are in GeV and GeV/c^2
@dataclass
class Particle:
    mat_type: MatterType
    particle: Hadron
    vector: FourVector
    mass: float

    # running `Particle | Particle` will shift the two particles into their mutual CoM frame
    def __or__(self, other):
        X = self.vector + other.vector
        V = X.get_velo()
        self.vector.shift_to_v(V)
        other.vector.shift_to_v(V)

    def __add__(self, other):
        return FourVector(self.vector.p + other.vector.p)

    def __sub__(self, other):
        return FourVector(self.vector.p - other.vector.p)


# A dataclass for managing events
# it has a unique ID, the list of involved particles, and the CoM frame
@dataclass
class ParticleEvent:
    id: int
    particles: np.ndarray


@dataclass
class Pairing:
    mat_type: MatterType
    particles: np.ndarray
    cm_velo: np.ndarray = field(init=False)
    pathway: int = field(default=0, init=False)

    def __post_init__(self):
        self.cm_velo = self.get_pair_v()

    def shift(self, V):
        self.particles[0].vector.shift_to_v(V)
        self.particles[1].vector.shift_to_v(V)

    def com(self):
        self.shift(self.cm_velo)

    def lab(self):
        self.shift(-self.cm_velo)

    def get_inter_particle(self):
        return self.particles[0] + self.particles[1]

    def get_pair_v(self):
        X = self.get_inter_particle()
        return X.get_velo()

    def set_pathway(self, pathway):
        self.pathway = pathway


@dataclass
class Deuteron(Particle):
    cm_velo: np.ndarray
    pathway: int
    particle: Hadron = field(default=Hadron.DEUTERON, init=False)
    mass: float = field(default=1.875612, init=False)

    def shift(self, V):
        self.vector.shift_to_v(V)

    def com(self):
        self.shift(self.cm_velo)

    def lab(self):
        self.shift(-self.cm_velo)


# A class to manage the complete catalogue of events
class EventCatalogue:
    def __init__(self, filename, debug=False):
        self.debug = debug
        self.events = []
        self.parse_data(filename)

    def parse_data(self, filename):
        with open(filename, "r") as f:
            content = f.read()

        events = content.split("\n\n")

        if self.debug is True:
            top = 100
        else:
            top = -1

        id = 0
        for event in events[1:top]:
            part_list = []
            particles = event.split("\n")
            for particle in particles:
                raw_data = particle.split()
                data = [int(raw_data[0]), *map(float, raw_data[1:])]
                vect = FourVector([data[1], data[2], data[3], data[4]])
                p_info = self.particle_decode(data[0])
                part_list.append(Particle(MatterType(p_info[0]), Hadron(p_info[1]), vect, data[5]))

            id += 1
            pl = np.array(part_list)
            self.events.append(ParticleEvent(id, pl))

    def particle_decode(self, value):
        match value:
            case 2212:
                return ("m", "p")
            case -2212:
                return ("a", "p")
            case 2112:
                return ("m", "n")
            case -2112:
                return ("a", "n")
            case _:
                raise ValueError(f"{value} is an unrecognised particle label; "
                                 "only (+/-)2212 and (+/-)2112 are valid.")


# %%
if __name__ == "__main__":
    pythia_EC = EventCatalogue("deuteron.dat")
