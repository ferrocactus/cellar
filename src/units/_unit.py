import logging
from abc import ABC, abstractmethod


class Unit(ABC):
    def __init__(self, name='Root'):
        self.name = name

    @abstractmethod
    def get(self):
        pass
