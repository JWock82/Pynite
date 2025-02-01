class LoadCombo():
    """A class that stores all the information necessary to define a load combination.
    """

    def __init__(self, name, combo_tags=None, factors={}):
        """Initializes a new load combination.

        :param name: A unique name for the load combination.
        :type name: str
        :param combo_tags: A list of tags for the load combination. This is a list of any strings you would like to use to categorize your load combinations. It is useful for separating load combinations into strength, service, or overstrength combinations as often required by building codes. This parameter has no effect on the analysis, but it can be used to restrict analysis to only the load combinations with the tags you specify.
        :type combo_tags: list, optional
        :param factors: A dictionary of load case names (`keys`) followed by their load factors (`items`). For example, the load combination 1.2D+1.6L would be represented as follows: `{'D': 1.2, 'L': 1.6}`. Defaults to {}.
        :type factors: dict, optional
        """
        
        self.name = name             # A unique user-defined name for the load combination
        self.combo_tags = combo_tags   # Used to categorize the load combination (e.g. strength or serviceability)
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
