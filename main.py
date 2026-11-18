# %%
import data_parser as dp
import combinations as cmb
from post_processing import get_cross_section
import numpy as np
import matplotlib.pyplot as plt


# if __name__ == "__main__":
# %%
EC = dp.EventCatalogue("deuteron.dat", debug=False)


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
            print(f"{itr}% complete...")
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


# %%
plt.plot(cutoffs, sigmas)
plt.xlabel("cutoff, $k$  $(GeV)$")
plt.ylabel(r"cross-section, $\sigma$  $(mb)$")
plt.title("Cross-Section against Cutoff")

# %%
all_ds = np.concatenate(all_ds)
print(all_ds)
# %%
