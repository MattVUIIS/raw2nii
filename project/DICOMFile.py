class DICOMFile(object):
    def __repr__(self):
        return 'DICOMFile' + repr(self.__dict__)

class ImageFrame(object):
    def __repr__(self):
        return 'ImageFrame' + repr(self.__dict__)

class MRFrame(object):
    def __repr__(self):
        return 'MRFrame' + repr(self.__dict__)
