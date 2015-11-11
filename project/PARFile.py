class PARFile:
    def __init__(self):
        #lambda function is used as an object to set arbitrary attributes on
        self.gen_info = lambda: None
        self.fields = []
        self.slices = None
        self.problemreading = False
        self.ResToolsVersion = None

    def __repr__(self):
        s = []
        for key, val in self.__dict__.items():
            s.append('{0}="{1}"'.format(key, val))
        return "<PARFile {0}>".format(" ".join(s))
