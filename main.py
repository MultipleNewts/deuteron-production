# %%
import data_parser as dp
import combinations as cmb


# %%
EC = dp.EventCatalogue("deuteron.dat", debug=True)


# %%
vcs = []
for event in EC.events:
    vcs.append(cmb.get_valid_pairs(event))

# %%
print(vcs[3])
# %%
