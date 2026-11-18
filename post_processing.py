# %%
def get_cross_section(deuterons, nevents):
    """
        Computes the cross-section given the number of deuterons and number of events

        Parameters
        ----------
        deuterons : `array`
            a flat numpy array containing all detected deuterons
        nevents : `int`
            the number of events tested

        Returns
        -------
        sigma : `float`
            the cross-section, sigma, of the test
    """
    count = deuterons.size
    return (count/nevents)*100e-3
