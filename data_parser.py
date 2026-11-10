# %%
from dataclasses import dataclass
import numpy as np
from four_vectors import FourVector


# A particle dataclass for the datastructure
# the particle string starts with either "m-" or "a-" for matter and antimatter respectively
# then it is succeeded by "-p" or "-n" for proton or neutron respectively
# it also contains the corresponds four-vector and mass
# Note: energy and mass-energy are in GeV and GeV/c^2
@dataclass
class Particle:
    particle: str
    vector: FourVector
    mass: float


# A dataclass for managing events
# it has a unique ID, the list of involved particles, and the CoM frame
@dataclass
class ParticleEvent:
    id: int
    particles: np.ndarray


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
            top = 10
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
                name = self.particle_decode(data[0])
                part_list.append(Particle(name, vect, data[5]))

            id += 1
            pl = np.array(part_list)
            self.events.append(ParticleEvent(id, pl))

    def particle_decode(self, value):
        match value:
            case 2212:
                return "mp"
            case -2212:
                return "ap"
            case 2112:
                return "mn"
            case -2112:
                return "an"
            case _:
                raise ValueError(f"{value} is an unrecognised particle label; "
                                 "only (+/-)2212 and (+/-)2112 are valid.")


# %%
if __name__ == "__main__":
    pythia_EC = EventCatalogue("deuteron.dat")
