#__all__ = []




class AnalysisMethod:
    '''
    Basic class for analysis methods. Inherited by SingleMethod and ComparisonMethod.
    '''
    
    def __init__(self, short_name, long_name, description):
        self.short_name = short_name
        self.long_name = long_name
        self.description = description

    def __str__(self):
        return "%s (%s): %s" % (self.long_name, self.short_name, self.description)

    def print_members(self):
        members = sorted([attr for attr in dir(self) if not callable(getattr(self,attr)) and not attr.startswith("__")])
        for m in members:
            print "%s = %s" % (m, getattr(self, m))

    def message(self, msg):
        print "[%s] %s" % (self.short_name, msg)


class SingleConditionMethod(AnalysisMethod):
    '''
    Class to be inherited by analysis methods that determine essentiality in a single condition (e.g. Gumble, Binomial, HMM).
    '''

    def __init__(self, short_name, long_name, description, ctrldata, replicates="Sum", normalization=None, LOESS=False, NTerminus=0.0, CTerminus=0.0):
        AnalysisMethod.__init__(self, short_name, long_name, description)
        self.ctrldata = ctrldata
        self.replicates = replicates
        self.normalization = normalization
        self.LOESS = LOESS
        self.NTerminus = NTerminus
        self.CTerminus = CTerminus


class DualConditionMethod(AnalysisMethod):
    '''
    Class to be inherited by analysis methods that determine changes in essentiality between two conditions (e.g. Resampling, DEHMM).
    '''

    def __init__(self, short_name, long_name, description, ctrldata, expdata, normalization, replicates="Sum", LOESS=False, NTerminus=0.0, CTerminus=0.0):
        AnalysisMethod.__init__(self, short_name, long_name, description)
        self.ctrldata = ctrldata
        self.expdata = expdata
        self.normalization = normalization
        self.replicates = replicates
        self.LOESS = LOESS
        self.NTerminus = NTerminus
        self.CTerminus = CTerminus





