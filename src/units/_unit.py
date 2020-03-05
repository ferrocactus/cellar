from abc import ABC, abstractmethod


class Unit(ABC):
    def __init__(self, verbose=False, name='NoName', **kwargs):
        self.__verbose = verbose
        self.__name = name
        self.kwargs = kwargs

    @abstractmethod
    def get(self):
        pass

    def vprint(self, *args):
        if self.__verbose == True:
            print(self.__name + ":", *args)