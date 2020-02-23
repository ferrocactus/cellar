from abc import ABC, abstractmethod


class Unit(ABC):
    def __init__(self, verbose=False, **args):
        self._verbose = verbose
        self.args = args

    @abstractmethod
    def get(self):
        pass

    def vprint(self, *args):
        if self._verbose == True:
            print(*args)