# %%
import data_parser as dp
import combinations as cmb


# if __name__ == "__main__":
# %%
EC = dp.EventCatalogue("deuteron.dat", debug=False)


# %%
ds = []
nevents = 100
itr = 1
for event in EC.events[:nevents]:
    if itr % (nevents//10) == 0:
        print(f"{itr}% complete...")
    event_ds = cmb.form_deuterons(event, cutoff=0.5)
    if event_ds is not None:
        ds.append(event_ds)
    itr += 1

# %%
print(ds)
# %%
