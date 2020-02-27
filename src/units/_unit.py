from abc import ABC, abstractmethod


class Unit(ABC):
    def __init__(self, verbose=False, **kwargs):
        self._verbose = verbose
        self.kwargs = kwargs
        self.name = 'NoName'

    @abstractmethod
    def get(self):
        pass

    def vprint(self, *args):
        if self._verbose == True:
            print(self.name + ":", *args)