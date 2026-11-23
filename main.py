# %%
import data_parser as dp
import combinations as cmb
from post_processing import get_cross_section, get_trans_momentum
import numpy as np
import matplotlib.pyplot as plt


# if __name__ == "__main__":
# %%
EC = dp.EventCatalogue("deuteron.dat", debug=False)
print(len(EC.events))

# %%
nevents = 50
cutoffs = np.arange(30, step=2.5)/10
sigmas = []
all_ds = []

for cutoff in cutoffs:
    print(f"The current cutoff is: {cutoff}")
    ds = []
    itr = 1
    for event in EC.events[:nevents]:
        if itr % (nevents//10) == 0:
            print(f"{itr*100/nevents}% complete...")
        event_ds = cmb.form_deuterons(event, cutoff=cutoff)
        if event_ds is not None:
            ds.append(event_ds)
        itr += 1
    if len(ds) != 0:
        flat_ds = np.concatenate(ds)
        sigmas.append(get_cross_section(flat_ds, nevents))
        all_ds.append(flat_ds)
    else:
        sigmas.append(0)
flat_ds = np.concatenate(all_ds)

# %%
plt.plot(cutoffs, sigmas)
plt.xlabel("Cutoff, $k$  $(GeV)$")
plt.ylabel(r"Cross-Section, $\sigma$  $(mb)$")
plt.title("Cross-Section against Cutoff")
plt.show()

# %%
tps = []
for deut in all_ds:
    tps.append(get_trans_momentum(deut))
# %%
plt.hist(tps, bins=50)
plt.xlabel("Transverse Momentum, $p_T$  $(GeV)$")
plt.ylabel(r"# of Occurrences, $n$")
plt.title("Distribution of Transverse Momenta")

# %%
nevents = len(EC.events)
sigmas = []
all_ds = []

ds = []
itr = 1
for event in EC.events[:nevents]:
    if itr % (nevents//1000) == 0:
        print(f"{itr*100/nevents:.4f}% complete...")
    event_ds = cmb.form_dr_deuterons(event)
    if event_ds is not None:
        ds.append(event_ds)
    itr += 1
if len(ds) != 0:
    flat_ds = np.concatenate(ds)
    sigmas.append(get_cross_section(flat_ds, nevents))
    all_ds.append(flat_ds)
else:
    sigmas.append(0)
flat_ds = np.concatenate(all_ds)
print("DONE.")

# %%
tps = []
for deut in flat_ds:
    tps.append(get_trans_momentum(deut))

paths = [1, 2, 3, 4]
all_path_tps = []
for path in paths:
    path_tps = []
    for deut in flat_ds:
        if deut.pathway[0] == path:
            path_tps.append(get_trans_momentum(deut))
    all_path_tps.append(path_tps)

# %%
plt.hist(tps, bins=50)
plt.xlabel("Transverse Momentum, $p_T$  $(GeV)$")
plt.ylabel(r"# of Occurrences, $n$")
plt.title("Distribution of Transverse Momenta")
plt.show()
# %%
plt.hist(all_path_tps[2], bins=30, alpha=0.5, label=r"$N_{1}N_{2}\to d\pi$")
plt.hist(all_path_tps[3], bins=30, alpha=0.5, label=r"$pn\to d\pi^{+}\pi^{-}$")
plt.hist(all_path_tps[0], bins=30, alpha=0.5, label=r"$pn\to d\gamma$")
plt.hist(all_path_tps[1], bins=30, alpha=0.5, label=r"$N N\to d\pi^{0}\pi$")
plt.xlabel("Transverse Momentum, $p_T$  $(GeV)$")
plt.ylabel(r"# of Occurrences, $n$")
plt.title("Distribution of Transverse Momenta")
plt.legend()
plt.show()

# %%
