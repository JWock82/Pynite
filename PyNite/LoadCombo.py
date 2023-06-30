class LoadCombo():
    """A class that stores all the information necessary to define a load combination.
    """

    def __init__(self, name, combo_type=None, factors={}):
        """Initializes a new load combination.

        :param name: A unique name for the load combination.
        :type name: str
        :param combo_type: The type of load combination. This can be any string you would like to use to categorize your load combinations. It is useful for separating load combinations into strength, service, or overstrength combinations as often required by building codes. This parameter has no effect on the analysis. It's simply a tool to help you filter results later on. Defaults to `None`.
        :type combo_type: str, optional
        :param factors: A dictionary of load case names (`keys`) followed by their load factors (`items`). For example, the load combination 1.2D+1.6L would be represented as follows: `{'D': 1.2, 'L': 1.6}`. Defaults to {}.
        :type factors: dict, optional
        """
        
        self.name = name             # A unique user-defined name for the load combination
        self.combo_type = combo_type # Used to track what type of load combination (e.g. strength or serviceability)
                                     # The 'combo_type' is just a placeholder for future features for now
        self.factors = factors       # A dictionary containing each load case name and associated load factor
    
    def AddLoadCase(self, case_name, factor):
        '''
        Adds a load case with its associated load factor
        '''

        self.factors[case_name] = factor
    
    def DeleteLoadCase(self, case_name):
        '''
        Deletes a load case with its associated load factor
        '''

        del self.factors[case_name]
