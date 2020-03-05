from abc import ABC, abstractmethod


class Unit(ABC):
    def __init__(self, verbose=False, name='NoName', **kwargs):
        self._verbose = verbose
        self.name = name
        self.kwargs = kwargs

    @abstractmethod
    def get(self):
        pass

    def vprint(self, *args):
        if self._verbose == True:
            print(self.name + ":", *args)