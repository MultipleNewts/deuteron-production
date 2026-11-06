# %%
from four_vectors import FourVector
from dataclasses import dataclass
import numpy as np


@dataclass
class Particle:
    particle: str
    vector: FourVector
    mass: float


@dataclass
class ParticleEvent:
    id: int
    particles: np.ndarray
    cm_frame: FourVector


class EventCatalogue:
    def __init__(self, filename):
        self.events = []
        self.parse_data(filename)

    def parse_data(self, filename):
        with open(filename, "r") as f:
            content = f.read()

        events = content.split("\n\n")

        id = 0
        for event in events[1:-1]:
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
            cm = FourVector([0, 0, 0, 0])  # placeholder value
            self.events.append(ParticleEvent(id, pl, cm))

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
pythia_EC = EventCatalogue("deuteron.dat")

# %%
print(pythia_EC.events[0].particles[0].type)
# %%
